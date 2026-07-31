#!/usr/bin/env python3
"""Create one evidence-ranked row per sample.

The report is descriptive, not diagnostic.  Small variants retain their
IthaGenes/ClinVar evidence.  Coverage-defined HBA gains are deliberately kept as
*structural candidates*: interval overlap cannot establish copy order, junction
sequence, or a named triplication allele.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path


def zygosity_from_gt(gt: str) -> str:
    alleles = (gt or "").strip().replace("|", "/").split("/")
    if alleles == ["1", "1"]:
        return "hom"
    if alleles in (["0", "1"], ["1", "0"]):
        return "het"
    if alleles == ["0", "0"]:
        return "ref"
    return ""


def tier(clinvar: str, functionality: str) -> int:
    classification = (clinvar or "").strip().lower()
    function = (functionality or "").strip().lower()
    if function == "causative" or "pathogenic" in classification:
        if classification.startswith("conflicting"):
            return 2
        return 1
    if classification.startswith("conflicting"):
        return 2
    if "uncertain" in classification or classification in (
        "",
        "other",
        "?",
        "no classification for the single variant",
    ):
        return 3
    if "benign" in classification:
        return 4
    return 3


TRANSCRIPT_TO_GENE = {
    "ENST00000335295": "HBB",
    "ENST00000251595": "HBA2",
    "ENST00000252242": "HBA1",
    "ENST00000380315": "HBB",
    "ENST00000866237": "HBA2",
    "ENST00000320868": "HBD",
    "ENST00001097508": "HBA1",
    "ENST00000485743": "HBB",
}


def to_hgvs(hgvs: str) -> str:
    """Convert transcript-prefixed c.HGVS to a gene-prefixed catalogue key."""

    value = str(hgvs)
    match = re.match(r"(ENST\d+)\.\d+:(c\..+)", value)
    if not match:
        return value
    gene = TRANSCRIPT_TO_GENE.get(match.group(1), match.group(1))
    return f"{gene}:{match.group(2)}"


def load_catalogue(path: str | Path) -> dict[str, tuple[str, str, str, str]]:
    """Return HGVS -> (classification, common name, functionality, source)."""

    catalogue: dict[str, tuple[str, str, str, str]] = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            hgvs = (row.get("HVGS") or "").strip()
            if not re.match(r"^[A-Z0-9]+:c\.", hgvs):
                continue
            source = (row.get("Source") or row.get("source") or "IthaGenes/ClinVar/HbVar").strip()
            catalogue[hgvs] = (
                (row.get("ClinVar classification") or "").strip(),
                (row.get("Common name") or "").strip(),
                (row.get("Functionality") or "").strip(),
                source,
            )
    return catalogue


def _excel_safe(value: object) -> str:
    text = "" if value is None else str(value)
    return "'" + text if text[:1] in ("-", "=", "+", "@") else text


def _display(common: str, hgvs: str) -> str:
    if common and hgvs:
        return f"{common} ({hgvs})"
    return common or hgvs


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) not in (3, 4):
        print(
            "usage: sample_report.py comprehensive_report.csv variants.csv "
            "output.csv [cnv_calls.csv]",
            file=sys.stderr,
        )
        return 2
    report_csv, variants_csv, out_csv = arguments[:3]
    cnv_csv = arguments[3] if len(arguments) == 4 else None

    catalogue = load_catalogue(variants_csv)
    cnv_by_sample: dict[str, list[dict[str, str]]] = {}
    if cnv_csv and os.path.exists(cnv_csv):
        with open(cnv_csv, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                cnv_by_sample.setdefault(row["Sample"], []).append(row)

    # Entries: tier, common_name, hgvs/evidence, zygosity, class, source.
    samples: dict[str, list[tuple[int, str, str, str, str, str]]] = {}
    with open(report_csv, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            sample = row["Sample"]
            samples.setdefault(sample, [])
            tool = row["Tool"]
            if tool == "Coverage":
                continue
            hgvs_raw = row.get("HGVS", "")
            consequence = row.get("Consequence", "")
            zygosity = zygosity_from_gt(row.get("Genotype", ""))

            if "SV_" in consequence or tool in ("Sniffles", "CuteSV"):
                name = re.split(r"\s*\(", str(hgvs_raw).strip())[0] if hgvs_raw else consequence
                if name.startswith("unknown"):
                    entry = (3, name, "", zygosity, "Unclassified SV", "caller evidence")
                else:
                    entry = (1, name, "", zygosity, "Causative", "IthaCNVs + read junction")
                if not any(existing[1] == entry[1] for existing in samples[sample]):
                    samples[sample].append(entry)
                continue

            hgvs = to_hgvs(hgvs_raw)
            clinvar, common, functionality, source = catalogue.get(
                hgvs, ("", "", "", "unmatched")
            )
            rank = tier(clinvar, functionality)
            classification = (clinvar or "").strip()
            if classification and classification not in (
                "?",
                "other",
                "not provided",
                "no classification for the single variant",
            ):
                klass = classification
            elif functionality.strip().lower() == "causative":
                klass = "Causative"
            else:
                klass = "unclassified"
            samples[sample].append((rank, common, hgvs, zygosity, klass, source))

    # Coverage-defined gains are candidates, not named causal alleles.
    for sample, rows in cnv_by_sample.items():
        samples.setdefault(sample, [])
        for row in rows:
            if row.get("Type") != "gain":
                continue
            candidate_name = (row.get("Name") or "").strip()
            naming_status = (row.get("Naming_status") or "").strip()
            details = f"{row.get('Region', '')}, ~{row.get('Fold', '?')}x depth"
            if candidate_name:
                details += f"; catalogue candidate {candidate_name}"
            common = "HBA copy-number gain candidate"
            klass = "Unresolved structural gain"
            source = "coverage"
            if naming_status:
                source += f"; {naming_status}"
            entry = (3, common, details, "", klass, source)
            if not any(existing[1] == common and existing[2] == details for existing in samples[sample]):
                samples[sample].append(entry)

    with open(out_csv, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Sample",
                "Result",
                "Common_name",
                "Primary_HGVS",
                "Primary_zygosity",
                "N_pathogenic",
                "N_conflicting",
                "N_vus",
                "N_benign",
                "All_variants_ranked",
            ]
        )
        for sample in sorted(samples):
            variants = sorted(samples[sample], key=lambda item: (item[0], item[1], item[2]))
            counts = {1: 0, 2: 0, 3: 0, 4: 0}
            for rank, *_ in variants:
                counts[rank] += 1

            primaries = [variant for variant in variants if variant[0] <= 2]
            primary_names = "; ".join(variant[1] for variant in primaries) if primaries else ""
            primary_hgvs = "; ".join(variant[2] for variant in primaries if variant[2])
            primary_zygosity = "; ".join(variant[3] for variant in primaries if variant[3])

            if counts[1]:
                result = "Causative variant detected"
            elif counts[2]:
                result = "Conflicting classification — review"
            elif counts[3]:
                result = "Unresolved/VUS finding — review"
            elif counts[4]:
                result = "Benign only"
            else:
                result = "No filtered variant detected"

            details = (
                "; ".join(
                    f"{_display(variant[1], variant[2])} "
                    f"[{variant[4]}{', ' + variant[3] if variant[3] else ''}; {variant[5]}]"
                    for variant in variants
                )
                if variants
                else "none"
            )
            writer.writerow(
                [
                    sample,
                    result,
                    _excel_safe(primary_names) or "none",
                    _excel_safe(primary_hgvs),
                    primary_zygosity,
                    counts[1],
                    counts[2],
                    counts[3],
                    counts[4],
                    details,
                ]
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
