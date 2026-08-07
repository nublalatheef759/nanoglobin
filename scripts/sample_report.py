#!/usr/bin/env python3
"""
sample_report.py -- one row per sample, unifying SNVs and structural variants
(deletions AND copy-number gains), ranked by clinical significance and flagged
with classification, source, and zygosity. Designed for high-volume scanning:
the left columns give a triage flag and the primary finding; the full ranked
variant list is pushed right, consulted only on demand.

Nothing is hidden: benign variants sort to the bottom but stay visible, because a
classification benign in one population/submitter may matter in another (and
ClinVar 'Conflicting' variants are surfaced explicitly).

The Result column is DESCRIPTIVE, not diagnostic: it reports what the pipeline
found (causative variant detected / conflicting / VUS only / none), leaving the
clinical conclusion to the clinician.

Inputs:
  comprehensive_report.csv   per-variant calls (SNVs + SV deletions)
  variants.csv               SNV/indel catalogue (ClinVar/HbVar classification)
  [cnv_calls.csv]            optional copy-number calls from detect_cnv.py (gains)

Output columns:
  Sample, Result, Common_name, Primary_HGVS, Primary_zygosity,
  N_pathogenic, N_conflicting, N_vus, N_benign, All_variants_ranked
"""
import csv, re, sys, os


def zygosity_from_gt(gt):
    g = (gt or "").strip().replace("|", "/")
    a = g.split("/")
    if a == ["1", "1"]:
        return "hom"
    if a in (["0", "1"], ["1", "0"]):
        return "het"
    if a == ["0", "0"]:
        return "ref"
    return ""


def tier(clinvar, functionality):
    c = (clinvar or "").strip().lower()
    f = (functionality or "").strip().lower()
    if f == "causative" or "pathogenic" in c:
        if c.startswith("conflicting"):
            return 2
        return 1
    if c.startswith("conflicting"):
        return 2
    if "uncertain" in c or c in ("", "other", "?",
                                 "no classification for the single variant"):
        return 3
    if "benign" in c:
        return 4
    return 3


TRANSCRIPT_TO_GENE = {
    "ENST00000335295": "HBB", "ENST00000251595": "HBA2", "ENST00000252242": "HBA1",
    "ENST00000380315": "HBB", "ENST00000866237": "HBA2", "ENST00000320868": "HBD",
    "ENST00001097508": "HBA1", "ENST00000485743": "HBB",
}


def to_hgvs(hgvs):
    """ENST00000335295.4:c.20A>T -> HBB:c.20A>T (gene-prefixed HGVS)."""
    s = str(hgvs)
    m = re.match(r"(ENST\d+)\.\d+:(c\..+)", s)
    if m:
        gene = TRANSCRIPT_TO_GENE.get(m.group(1), m.group(1))
        return "%s:%s" % (gene, m.group(2))
    return s


def load_catalogue(path):
    """gene:c. HGVS -> (clinvar_class, common_name, functionality). Skips junk rows."""
    cat = {}
    with open(path, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            hvgs = (r.get("HVGS") or "").strip()
            if not re.match(r"^[A-Z0-9]+:c\.", hvgs):
                continue
            cat[hvgs] = (
                (r.get("ClinVar classification") or "").strip(),
                (r.get("Common name") or "").strip(),
                (r.get("Functionality") or "").strip(),
            )
    return cat



def _excel_safe(v):
    """Excel reads a leading -, =, + or @ as a formula (e.g. -α3.7 -> #NAME?).
    Prefix a leading apostrophe so Excel renders the literal text. Only applied
    to standalone name cells, not the combined detail column."""
    v = "" if v is None else str(v)
    return "'" + v if v[:1] in ("-", "=", "+", "@") else v


def main():
    report_csv = sys.argv[1]
    variants_csv = sys.argv[2]
    out_csv = sys.argv[3]
    cnv_csv = sys.argv[4] if len(sys.argv) > 4 else None

    cat = load_catalogue(variants_csv)

    cnv_by_sample = {}
    if cnv_csv and os.path.exists(cnv_csv):
        with open(cnv_csv, encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                cnv_by_sample.setdefault(r["Sample"], []).append(r)

    # each entry: (tier, common_name, hgvs, zyg, klass, source)
    samples = {}
    with open(report_csv, encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            s = row["Sample"]
            samples.setdefault(s, [])
            tool = row["Tool"]
            if tool == "Coverage":
                continue
            hgvs_raw = row.get("HGVS", "")
            cons = row.get("Consequence", "")
            zyg = zygosity_from_gt(row.get("Genotype", ""))

            if "SV_" in cons or tool in ("Sniffles", "CuteSV"):
                # structural deletion: name already the identify_sv result
                name = re.split(r"\s*\(", str(hgvs_raw).strip())[0] if hgvs_raw else cons
                if name.startswith("unknown"):
                    common, hgvs_disp, klass, src, t = name, "", "Unclassified SV", "-", 3
                else:
                    common, hgvs_disp, klass, src, t = name, "", "Causative", "IthaCNVs", 1
                # dedup: Sniffles + CuteSV report the same deletion
                if any(e[1] == common for e in samples[s]):
                    continue
                samples[s].append((t, common, hgvs_disp, zyg, klass, src))
            else:
                hgvs = to_hgvs(hgvs_raw)
                clinvar, common, func = cat.get(hgvs, ("", "", ""))
                t = tier(clinvar, func)
                cv = (clinvar or "").strip()
                if cv and cv not in ("?", "other", "not provided",
                                     "no classification for the single variant"):
                    klass = cv
                elif func.strip().lower() == "causative":
                    klass = "Causative"
                else:
                    klass = "unclassified"
                samples[s].append((t, common, hgvs, zyg, klass, "ClinVar/HbVar"))

    # copy-number GAINS from detect_cnv (losses already come via SV callers)
    for s_name, rows in cnv_by_sample.items():
        samples.setdefault(s_name, [])
        for r in rows:
            if r["Type"] == "gain":
                named = (r.get("Name") or "").strip()
                # Coverage establishes that extra material is present; it does NOT
                # establish copy order, chromosome assignment, hybrid orientation or
                # the junction. A reciprocal catalogue overlap is therefore reported
                # as a candidate match, never as a definitive named triplication.
                if named:
                    match = re.split(r"\s*\(\d+%", named)[0].strip()
                    common = "unresolved structural-gain candidate (best catalogue match: %s)" % named
                    klass = "gain, unresolved (coverage only)"
                else:
                    common = "unresolved structural-gain candidate"
                    klass = "gain, unresolved (coverage only)"
                hgvs_disp = "%s, ~%s× depth" % (r["Region"], r.get("Fold", "?"))
                if not any(e[1] == common for e in samples[s_name]):
                    samples[s_name].append((1, common, hgvs_disp, "", klass, "coverage"))

    def display(common, hgvs):
        """Combined 'Name (HGVS)' for the detail column."""
        if common and hgvs:
            return "%s (%s)" % (common, hgvs)
        return common or hgvs

    with open(out_csv, "w", newline="", encoding="utf-8-sig") as out:
        w = csv.writer(out)
        w.writerow(["Sample", "Result", "Common_name", "Primary_HGVS",
                    "Primary_zygosity", "N_pathogenic", "N_conflicting",
                    "N_vus", "N_benign", "All_variants_ranked"])
        for s in sorted(samples):
            variants = sorted(samples[s], key=lambda x: x[0])
            n = {1: 0, 2: 0, 3: 0, 4: 0}
            for t, *_ in variants:
                n[t] += 1

            # primary = tier 1/2 findings
            primaries = [v for v in variants if v[0] <= 2]
            prim_names = "; ".join(v[1] for v in primaries) if primaries else ""
            prim_hgvs = "; ".join(v[2] for v in primaries if v[2]) if primaries else ""
            prim_zyg = "; ".join(v[3] for v in primaries if v[3]) if primaries else ""

            # descriptive Result flag (NOT a diagnosis)
            if n[1]:
                result = "Causative variant detected"
            elif n[2]:
                result = "Conflicting classification — review"
            elif n[3]:
                result = "VUS only — review"
            elif n[4]:
                result = "Benign only"
            else:
                result = "No variant detected"

            all_v = "; ".join(
                "%s [%s%s]" % (display(v[1], v[2]), v[4], (", " + v[3]) if v[3] else "")
                for v in variants
            ) if variants else "none"

            w.writerow([s, result, _excel_safe(prim_names) or "none", _excel_safe(prim_hgvs), prim_zyg,
                        n[1], n[2], n[3], n[4], all_v])


if __name__ == "__main__":
    main()
    
    

    