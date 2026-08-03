"""Score final report rows against a declared synthetic truth manifest.

This utility is for software regression fixtures. It reads
``results/comprehensive_report.csv`` and matches reported SNV/SV rows to the truth
VCFs listed in ``simulation/benchmark_truth.tsv``. It does not provide biological
or clinical validation.

Usage::

    python scripts/benchmark.py \
      --report results/comprehensive_report.csv \
      --config simulation/benchmark_truth.tsv \
      --sv-tool Sniffles \
      --out results/synthetic_benchmark.tsv
"""
import argparse, csv, re, sys

RECIP, SIZE_TOL = 0.5, 0.20
SVLEN_RE = re.compile(r"SV_[A-Z]+_(-?\d+)bp")


def parse_truth(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            is_sv = "SVTYPE" in info
            if is_sv:
                start = int(f[1])
                end = int(info.get("END", start + abs(int(info.get("SVLEN", 0)))))
                out.append({"sv": True, "chrom": f[0], "start": start, "end": end})
            else:
                out.append({"sv": False, "chrom": f[0], "pos": int(f[1]),
                            "ref": f[3], "alt": f[4]})
    return out


def load_report(path, sv_tool):
    """Return {sample: {'sv': [...], 'snv': [...]}}."""
    samples = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            s = samples.setdefault(r["Sample"], {"sv": [], "snv": []})
            tool = r["Tool"]
            if tool in ("Sniffles", "CuteSV"):
                if tool != sv_tool:
                    continue
                m = SVLEN_RE.search(r.get("Consequence", ""))
                if not m:
                    continue
                start = int(r["Position"])
                size = abs(int(m.group(1)))
                s["sv"].append({"chrom": r["Chromosome"], "start": start,
                                "end": start + size})
            elif tool == "Clair3":
                s["snv"].append({"chrom": r["Chromosome"], "pos": int(r["Position"]),
                                 "ref": r["Ref"], "alt": r["Alt"]})
    return samples


def sv_match(t, c):
    if t["chrom"] != c["chrom"]:
        return False
    tl, cl = t["end"] - t["start"], c["end"] - c["start"]
    if tl <= 0 or cl <= 0:
        return False
    ov = max(0, min(t["end"], c["end"]) - max(t["start"], c["start"]))
    if min(ov / tl, ov / cl) < RECIP:
        return False
    return (1 - SIZE_TOL) <= (cl / tl) <= (1 + SIZE_TOL)


def snv_match(t, c):
    return (t["chrom"] == c["chrom"] and t["pos"] == c["pos"]
            and t["ref"] == c["ref"] and t["alt"] == c["alt"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--sv-tool", default="Sniffles")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    report = load_report(a.report, a.sv_tool)
    rows = []
    with open(a.config) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sample, cls = r["sample"], r["class"]
            calls = report.get(sample, {"sv": [], "snv": []})
            pool = calls["sv"] if cls == "SV" else calls["snv"]
            matcher = sv_match if cls == "SV" else snv_match

            if r["truth_vcf"] == "-":
                tp, fn, fp, detected = 0, 0, len(pool), "-"
            else:
                truths = [t for t in parse_truth(r["truth_vcf"])
                          if t["sv"] == (cls == "SV")]
                tp = fn = 0
                used = set()
                for t in truths:
                    hit = next((i for i, c in enumerate(pool)
                                if i not in used and matcher(t, c)), None)
                    if hit is not None:
                        tp += 1; used.add(hit)
                    else:
                        fn += 1
                fp = len(pool) - len(used)
                detected = "yes" if tp else "no"

            prec = tp / (tp + fp) if (tp + fp) else None
            rec = tp / (tp + fn) if (tp + fn) else None
            rows.append({"sample": sample, "class": cls, "TP": tp, "FP": fp, "FN": fn,
                         "precision": "%.3f" % prec if prec is not None else "NA",
                         "recall": "%.3f" % rec if rec is not None else "NA",
                         "detected": detected})

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample","class","TP","FP","FN","precision","recall","detected"], delimiter="\t")
        w.writeheader(); w.writerows(rows)

    cols = ["sample","class","TP","FP","FN","precision","recall","detected"]
    wd = {c: max(len(c), max((len(str(r[c])) for r in rows), default=0)) for c in cols}
    print(" | ".join(c.ljust(wd[c]) for c in cols))
    print("-+-".join("-"*wd[c] for c in cols))
    for r in rows:
        print(" | ".join(str(r[c]).ljust(wd[c]) for c in cols))


if __name__ == "__main__":
    main()
    
    
