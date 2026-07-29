#!/usr/bin/env python3
"""
Build an IthaGenes-primary naming layer, cross-referenced with ClinVar
pathogenicity and gnomAD frequency.

Design (Option A):
  SPINE   = IthaGenes (curated, haemoglobinopathy-specific, causative-enriched,
            carries traditional Hb names + functionality + phenotype)
  JOIN 1  = ClinVar classification (from existing variants.csv) on HGVS
  JOIN 2  = gnomAD allele frequencies (from existing variants.csv) on HGVS
  JOIN 3  = cohort frequency (from existing variants.csv) on HGVS

Output: naming_layer.csv  — every IthaGenes variant, enriched where a ClinVar/
gnomAD/cohort match exists. Nothing from variants.csv is lost that isn't in
IthaGenes: those rows are appended at the end tagged source='ClinVar-only' so
the layer is a strict superset (no coverage regression).

This is a STANDALONE artifact. It does not touch the pipeline or any validated
output. Switching the pipeline to use it is a separate, later step.
"""
import argparse
import sys
import pandas as pd
import numpy as np


def norm_hgvs(s):
    """Normalise an HGVS string for joining: strip, collapse internal spaces,
    split pipe-joined multi-gene rows are handled upstream (explode)."""
    if pd.isna(s):
        return ""
    s = str(s).strip()
    # collapse whitespace around ':' and internal double spaces
    s = " ".join(s.split())
    s = s.replace(": c.", ":c.").replace(": g.", ":g.").replace(": p.", ":p.")
    return s


def explode_pipes(df, col):
    """IthaGenes lists multi-gene variants as 'HBA1:c.X | HBA2:c.X' in one row.
    Explode into one row per HGVS so joins work, keeping all other columns."""
    df = df.copy()
    df[col] = df[col].astype(str).str.split(r"\s*\|\s*")
    df = df.explode(col)
    df[col] = df[col].map(norm_hgvs)
    df = df[df[col] != ""]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ithagenes", default="IthaGenes_full.csv")
    ap.add_argument("--variants", default="variants.csv")
    ap.add_argument("--out", default="naming_layer.csv")
    ap.add_argument("--genes", default="HBA,HBB",
                    help="comma-separated gene prefixes to keep (default HBA,HBB)")
    args = ap.parse_args()

    keep_prefixes = tuple(g.strip() for g in args.genes.split(","))

    # ---- load IthaGenes (the spine) ----
    itha = pd.read_csv(args.ithagenes, dtype=str, keep_default_na=False)
    itha.columns = [c.strip() for c in itha.columns]
    # standardise the HGVS column name
    hgvs_col = [c for c in itha.columns if c.lower().replace(" ", "") in
                ("hgvsname", "hgvs")][0]
    itha = itha.rename(columns={hgvs_col: "HGVS"})
    itha = explode_pipes(itha, "HGVS")

    # keep only HBA/HBB (the genes the pipeline calls)
    mask = itha["HGVS"].str.startswith(keep_prefixes)
    itha = itha[mask].copy()
    itha = itha.drop_duplicates(subset="HGVS")
    print(f"IthaGenes HBA/HBB variants (spine): {len(itha)}", file=sys.stderr)

    # ---- load existing variants.csv (ClinVar + gnomAD + cohort) ----
    var = pd.read_csv(args.variants, dtype=str, keep_default_na=False)
    var.columns = [c.strip() for c in var.columns]
    # the HGVS column in variants.csv is 'HVGS' (note the typo in the original)
    vhgvs = [c for c in var.columns if c.upper() in ("HVGS", "HGVS")][0]
    var = var.rename(columns={vhgvs: "HGVS"})
    var = explode_pipes(var, "HGVS")
    var = var.drop_duplicates(subset="HGVS")

    # pick the cross-reference columns we want to carry over (only those present)
    xref_candidates = {
        "ClinVar classification": "clinvar_classification",
        "clinvar_review_status": "clinvar_review_status",
        "clinvar_id": "clinvar_id",
        "In gnomAD? (Yes/No)": "in_gnomad",
        "gnomAD total AF (%)": "gnomad_af_total",
        "gnomAD ME AF (%)": "gnomad_af_me",
        "Cohort_n": "cohort_n",
        "Cohort frequency (%)": "cohort_freq_pct",
        "chrom": "chrom", "pos": "pos", "ref": "ref", "alt": "alt",
    }
    have = {k: v for k, v in xref_candidates.items() if k in var.columns}
    xref = var[["HGVS"] + list(have.keys())].rename(columns=have)

    # ---- JOIN: IthaGenes spine LEFT JOIN cross-references ----
    merged = itha.merge(xref, on="HGVS", how="left")

    # tag provenance
    merged["in_ithagenes"] = "Yes"
    merged["in_clinvar"] = np.where(
        merged.get("clinvar_classification", pd.Series([""] * len(merged))).fillna("") != "",
        "Yes", "No")

    # ---- append ClinVar-only variants (present in variants.csv, absent from IthaGenes) ----
    # so the layer is a strict superset — no coverage regression
    itha_hgvs = set(itha["HGVS"])
    var_hba_hbb = var[var["HGVS"].str.startswith(keep_prefixes)].copy()
    clinvar_only = var_hba_hbb[~var_hba_hbb["HGVS"].isin(itha_hgvs)].copy()

    # build clinvar-only rows in the same schema
    co = pd.DataFrame({"HGVS": clinvar_only["HGVS"]})
    # carry the same xref columns
    for src, dst in have.items():
        co[dst] = clinvar_only[src].values
    # fill IthaGenes-native columns as blank
    for c in itha.columns:
        if c != "HGVS" and c not in co.columns:
            co[c] = ""
    co["in_ithagenes"] = "No"
    co["in_clinvar"] = np.where(
        co.get("clinvar_classification", pd.Series([""] * len(co))).fillna("") != "",
        "Yes", "No")

    # align columns and concat
    final = pd.concat([merged, co], ignore_index=True, sort=False)

    # order columns sensibly: identity first, then annotations
    front = ["HGVS"]
    for c in ["Common Name", "Hb Name", "Genes", "Functionality", "Phenotype",
              "Locus", "Position", "IthaID"]:
        if c in final.columns:
            front.append(c)
    rest = [c for c in final.columns if c not in front]
    final = final[front + rest]

    final.to_csv(args.out, index=False)

    # ---- report ----
    n_total = len(final)
    n_itha = (final["in_ithagenes"] == "Yes").sum()
    n_clinvar_only = (final["in_ithagenes"] == "No").sum()
    n_both = ((final["in_ithagenes"] == "Yes") & (final["in_clinvar"] == "Yes")).sum()
    n_itha_no_clinvar = ((final["in_ithagenes"] == "Yes") & (final["in_clinvar"] == "No")).sum()
    print("=" * 50, file=sys.stderr)
    print(f"naming_layer.csv written: {n_total} variants", file=sys.stderr)
    print(f"  IthaGenes-sourced (spine):        {n_itha}", file=sys.stderr)
    print(f"    of which also have ClinVar:     {n_both}", file=sys.stderr)
    print(f"    IthaGenes-only (no ClinVar):    {n_itha_no_clinvar}", file=sys.stderr)
    print(f"  ClinVar-only (not in IthaGenes):  {n_clinvar_only}", file=sys.stderr)
    print("=" * 50, file=sys.stderr)
    # functionality breakdown (the clinical value-add)
    if "Functionality" in final.columns:
        print("Functionality (IthaGenes-sourced):", file=sys.stderr)
        fb = final.loc[final["in_ithagenes"] == "Yes", "Functionality"].replace("", "(blank)").value_counts()
        for k, v in fb.items():
            print(f"    {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
