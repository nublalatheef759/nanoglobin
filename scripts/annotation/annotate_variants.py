"""
nano-thal: Clinical annotation of thalassaemia variants
Matches pipeline output against curated variant database
and generates clinical interpretation report.
"""

import pandas as pd
import sys
import os
import re
import csv

# β⁰/β⁺ severity classification, loaded once from databases/beta_classification.csv
_BETA_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "databases", "beta_classification.csv")

_BETA = {}   # normalised hgvs -> label
def _norm_beta(c):
    return re.sub(r'(del|dup|ins)[ACGT]+$', r'\1', c) if c else c
try:
    with open(_BETA_PATH, encoding="utf-8-sig") as _f:
        for _r in csv.DictReader(_f):
            _h = _norm_beta((_r.get("hgvs") or "").strip())
            if _h:
                _BETA[_h] = (_r.get("label") or "").strip()
except Exception:
    pass
    
from datetime import datetime

# === CONFIGURATION ===
TRANSCRIPT_TO_GENE = {
    "ENST00000335295.4": "HBB",
    "ENST00000251595.11": "HBA2",
    "ENST00000252242.10": "HBA1"
}

# === FUNCTIONS ===

def convert_hgvs(hgvs):
    """Convert transcript HGVS to gene HGVS (e.g., ENST00000335295.4:c.92+5G>C -> HBB:c.92+5G>C)"""
    if pd.isna(hgvs) or hgvs == "unknown":
        return None
    for transcript, gene in TRANSCRIPT_TO_GENE.items():
        if transcript in hgvs:
            return hgvs.replace(transcript, gene)
    return hgvs

def load_database(db_path):
    """Load the variant lookup database"""
    db = pd.read_csv(db_path)
    # Clean column names
    db.columns = [c.strip() for c in db.columns]
    return db

def annotate_snv_indels(variants_df, db):
    """Match SNV/indel variants against the database"""
    variants_df = variants_df.copy()
    
    # Convert HGVS format
    variants_df['HGVS_gene'] = variants_df['HGVS'].apply(convert_hgvs)
    
    # Extract just the c. notation for matching
    variants_df['c_notation'] = variants_df['HGVS_gene'].apply(
        lambda x: x.split(':')[1] if pd.notna(x) and isinstance(x, str) and ':' in x else None
    )
    
    # Try matching on HGVS
    db['c_notation'] = db['HVGS'].apply(
        lambda x: x.split(':')[1] if pd.notna(x) and isinstance(x, str) and ':' in x else str(x) if pd.notna(x) else None
    )
    
    annotated = variants_df.merge(
        db[['HVGS', 'Common name', 'Gene', 'Cohort_n', 'Cohort frequency (%)',
            'Functionality', 'In Ithagenes? (Yes/No)', 'In gnomAD? (Yes/No)',
            'gnomAD total AF (%)', 'gnomAD ME AF (%)']],
        left_on='HGVS_gene',
        right_on='HVGS',
        how='left',
        suffixes=('_pipeline', '_db')
    )
    
    return annotated

def _smart_split(cell):
    """Split a patient_summary cell into individual variants on '; ', but NOT on
    the '; ' that appears inside an identify_sv annotation like
    '(100% reciprocal; ithaID=2231)'. Paren-depth aware."""
    parts, depth, buf, i, s = [], 0, "", 0, str(cell)
    while i < len(s):
        c = s[i]
        if c == '(':
            depth += 1; buf += c
        elif c == ')':
            depth -= 1; buf += c
        elif s[i:i+2] == '; ' and depth == 0:
            parts.append(buf); buf = ""; i += 2; continue
        else:
            buf += c
        i += 1
    if buf:
        parts.append(buf)
    return parts


def annotate_structural(patient_summary_df, cnv_db):
    """Match structural variants against the IthaCNVs catalogue (cnvs.csv).

    Handles BOTH HBA and HBB structural variants. Each variant name in
    patient_summary is already the identify_sv match, e.g.
    '-α3.7 (type III) (100% reciprocal; ithaID=2231)'. We look each up by its
    ithaID (parsed from that name) against cnvs.csv, which is coordinate-derived,
    gene-agnostic, and includes HBB entries (Corfu, Lepore, IVS deletions).

    Replaces the former hardcoded sv_lookup dict (HBA-only, 8 entries) with the
    258-CNV catalogue. The HBB path is exercised by cnvs.csv but untested in the
    current dataset (no HBB SVs in the test samples).
    """
    import re
    results = []

    def lookup(var):
        m = re.search(r'ithaID=(\w+)', str(var))
        if m is not None:
            hit = cnv_db[cnv_db['ithaID'].astype(str) == m.group(1)]
            if not hit.empty:
                return hit.iloc[0]
        base = re.split(r'\s*\(', str(var).strip())[0].strip()
        hit = cnv_db[cnv_db['name'].astype(str).str.strip() == base]
        if not hit.empty:
            return hit.iloc[0]
        return None

    def emit(sample, gene, var, zyg):
        match = lookup(var)
        if match is not None:
            results.append({
                'Sample': sample, 'Variant_pipeline': var,
                'Variant_db': match.get('name', var),
                'Common_name': match.get('name', ''),
                'Gene': gene, 'Zygosity': zyg, 'Type': 'Structural',
                'Functionality': match.get('functionality', ''),
                'Locus': match.get('locus', ''),
                'HGVS': match.get('hgvs', ''),
                'IthaID': match.get('ithaID', ''),
                'Source': match.get('source', ''),
            })
        else:
            results.append({
                'Sample': sample, 'Variant_pipeline': var, 'Variant_db': var,
                'Common_name': 'NOT IN CATALOGUE', 'Gene': gene, 'Zygosity': zyg,
                'Type': 'Structural', 'Functionality': 'Unknown',
                'Locus': '', 'HGVS': '', 'IthaID': '', 'Source': '',
            })

    for _, row in patient_summary_df.iterrows():
        sample = row['Sample']
        for gene in ('HBA', 'HBB'):
            cell = row.get(gene)
            if not pd.notna(cell) or str(cell).strip() in ('', 'wild-type'):
                continue
            variants = _smart_split(cell)
            zygs = _smart_split(row.get(gene + '_zygosity', '') or '')
            for i, var in enumerate(variants):
                var = var.strip()
                if not var or var == 'wild-type':
                    continue
                if var.startswith('ENST') or ':c.' in var:
                    continue
                zyg = zygs[i].strip() if i < len(zygs) else 'unknown'
                emit(sample, gene, var, zyg)

    return pd.DataFrame(results)

def _is_causative(variant_name, gene, db):
    """True if this variant is causative/pathogenic per the database."""
    if not variant_name or str(variant_name).strip().lower() in ("", "wild-type", "nan"):
        return False
    base = re.split(r"\s*\(", str(variant_name).strip())[0].strip()
    for _, r in db.iterrows():
        if str(r.get("Gene", "")).upper() not in gene.upper():
            continue
        hvgs = str(r.get("HVGS", "")).strip()
        common = str(r.get("Common name", "")).strip()
        if base and (base == hvgs or base == common or base in hvgs):
            func = str(r.get("Functionality", "")).strip().lower()
            clin = str(r.get("ClinVar classification", "")).strip().lower()
            if func.startswith("non") or "benign" in clin:
                return False
            return func == "causative" or "pathogenic" in clin
    return False


def _alpha_dose(zygosity):
    z = str(zygosity).strip().lower()
    if "hom" in z or "comp" in z:
        return "hom"
    if "het" in z:
        return "het"
    return "none"


def check_coinheritance(patient_summary_df, db):
    """Flag patients co-inheriting CAUSATIVE HBA and HBB variants.

    Co-inherited alpha-thalassaemia suppresses HbA2, the diagnostic marker for
    beta-thal trait, so an HBB carrier can screen normal on HPLC. Benign variants
    do not do this, so both sides must be causative before flagging.
    """
    flags = []
    for _, row in patient_summary_df.iterrows():
        hba, hbb = row.get("HBA", ""), row.get("HBB", "")
        if not _is_causative(hba, "HBA", db) or not _is_causative(hbb, "HBB", db):
            continue
        dose = _alpha_dose(row.get("HBA_zygosity", ""))
        pct = COHORT_BELOW_CUTOFF[dose]
        flags.append({
            "Sample": row["Sample"],
            "HBA": hba,
            "HBA_zygosity": row.get("HBA_zygosity", ""),
            "HBB": hbb,
            "HBB_zygosity": row.get("HBB_zygosity", ""),
            "Flag": (
                "CO-INHERITED CAUSATIVE HBA + HBB VARIANTS. Co-inherited "
                "alpha-thalassaemia suppresses HbA2, the marker for beta-thal "
                "trait. In this cohort %.1f%% of HBB carriers with this HBA "
                "status fell below the 3.5%% HbA2 cutoff (vs 25.5%% with no HBA "
                "variant). HPLC may report normal; molecular confirmation indicated."
                % pct
            ),
            "Cohort_pct_below_cutoff": pct,
            "Risk_level": "HIGH" if dose == "hom" else "MODERATE",
        })
    return pd.DataFrame(flags)
    
def classify_mutation_type(hgvs, gene):
    """Classify HBB variants by mutation type, from databases/beta_classification.csv."""
    if gene != 'HBB':
        return 'N/A'
    if pd.isna(hgvs):
        return 'Unknown'
    c_part = hgvs.split(':')[1] if ':' in str(hgvs) else str(hgvs)
    return _BETA.get(_norm_beta(c_part), 'Unclassified')

def generate_report(sample, snv_results, sv_results, coinheritance_flags):
    """Generate a clinical report for a single sample"""
    report = []
    report.append(f"{'='*60}")
    report.append(f"NANO-THAL CLINICAL REPORT")
    report.append(f"{'='*60}")
    report.append(f"Sample: {sample}")
    report.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append(f"{'='*60}")
    
    # Structural variants
    sample_sv = sv_results[sv_results['Sample'] == sample]
    if not sample_sv.empty:
        report.append(f"\nSTRUCTURAL VARIANTS:")
        report.append(f"{'-'*40}")
        for _, row in sample_sv.iterrows():
            report.append(f"  Variant: {row['Variant_db']}")
            report.append(f"  Common name: {row['Common_name']}")
            report.append(f"  Gene: {row['Gene']}")
            report.append(f"  Zygosity: {row['Zygosity']}")
            report.append(f"  Functionality: {row['Functionality']}")
            if str(row.get('Locus','')).strip():
                report.append(f"  Locus: {row['Locus']}")
            if str(row.get('HGVS','')).strip():
                report.append(f"  HGVS: {row['HGVS']}")
            if str(row.get('IthaID','')).strip():
                report.append(f"  IthaID: {row['IthaID']}")
            report.append("")
    
    # SNV/Indel variants
    sample_snv = snv_results[snv_results['Sample'] == sample]
    if not sample_snv.empty:
        report.append(f"\nSNV/INDEL VARIANTS:")
        report.append(f"{'-'*40}")
        for _, row in sample_snv.iterrows():
            report.append(f"  Position: {row.get('Chromosome', 'N/A')}:{row.get('Position', 'N/A')}")
            hgvs_val = row.get('HGVS_gene', row.get('HGVS', 'N/A'))
            report.append(f"  HGVS: {hgvs_val if pd.notna(hgvs_val) else 'N/A'}")
            if pd.notna(row.get('Common name')):
                report.append(f"  Common name: {row['Common name']}")
            report.append(f"  Gene: {row.get('Gene_pipeline', row.get('Gene', 'N/A'))}")
            report.append(f"  Consequence: {row.get('Consequence', 'N/A')}")
            report.append(f"  Genotype: {row.get('Genotype', 'N/A')}")
            report.append(f"  Quality: {row.get('Quality', 'N/A')}")
            
            mutation_type = classify_mutation_type(
                row.get('HGVS_gene', ''), 
                row.get('Gene_pipeline', row.get('Gene', ''))
            )
            report.append(f"  Mutation type: {mutation_type}")
            
            if pd.notna(row.get('Functionality')):
                report.append(f"  Functionality: {row['Functionality']}")
            if pd.notna(row.get('In gnomAD? (Yes/No)')):
                report.append(f"  In gnomAD: {row['In gnomAD? (Yes/No)']}")
            report.append("")
    
    # Co-inheritance flags
    if not coinheritance_flags.empty:
            sample_flags = coinheritance_flags[coinheritance_flags['Sample'] == sample]
    else:
        sample_flags = coinheritance_flags
    if not sample_flags.empty:
        report.append(f"\n{'!'*60}")
        report.append(f"⚠ CO-INHERITANCE ALERT")
        report.append(f"{'!'*60}")
        for _, flag in sample_flags.iterrows():
            report.append(f"  {flag['Flag']}")
            report.append(f"  Risk level: {flag['Risk_level']}")
            report.append(f"  HBA: {flag['HBA']} ({flag['HBA_zygosity']})")
            report.append(f"  HBB: {flag['HBB']} ({flag['HBB_zygosity']})")
    
    report.append(f"\n{'='*60}")
    report.append(f"END OF REPORT")
    report.append(f"{'='*60}")
    
    return '\n'.join(report)


# === MAIN ===
if __name__ == "__main__":
    # Paths
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    variants_path = os.path.join(base_dir, "results", "all_variants_annotated.csv")
    patient_path = os.path.join(base_dir, "results", "patient_summary.csv")
    db_path = os.path.join(base_dir, "databases", "variants.csv")
    output_dir = os.path.join(base_dir, "results", "clinical_reports")
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("nano-thal Clinical Annotation Pipeline")
    print("=" * 40)
    
    # Load data
    print("Loading variant data...")
    variants = pd.read_csv(variants_path)
    patients = pd.read_csv(patient_path)
    db = load_database(db_path)
    cnv_db = load_database(os.path.join(base_dir, "databases", "cnvs.csv"))
    
    print(f"  Samples: {len(patients)}")
    print(f"  SNV/Indel variants: {len(variants)}")
    print(f"  Database variants: {len(db)}")
    
    # Annotate
    print("\nAnnotating SNV/indels...")
    snv_annotated = annotate_snv_indels(variants, db)
    
    print("Annotating structural variants...")
    sv_annotated = annotate_structural(patients, cnv_db)
    
    print("Checking co-inheritance...")
    coinheritance = check_coinheritance(patients, db)
    
    # Save annotated results
    snv_annotated.to_csv(os.path.join(output_dir, "snv_annotated.csv"), index=False)
    sv_annotated.to_csv(os.path.join(output_dir, "sv_annotated.csv"), index=False)
    coinheritance.to_csv(os.path.join(output_dir, "coinheritance_flags.csv"), index=False)
    
    # Generate per-sample reports
    print("\nGenerating clinical reports...")
    all_samples = patients['Sample'].unique()
    
    for sample in all_samples:
        report = generate_report(sample, snv_annotated, sv_annotated, coinheritance)
        
        report_path = os.path.join(output_dir, f"{sample}_report.txt")
        with open(report_path, 'w') as f:
            f.write(report)
        print(f"  Generated: {sample}_report.txt")
    
    # Summary
    print(f"\n{'='*40}")
    print(f"SUMMARY")
    print(f"{'='*40}")
    print(f"Reports generated: {len(all_samples)}")
    print(f"Co-inheritance alerts: {len(coinheritance)}")
    print(f"Database matches (SNV): {snv_annotated['HVGS'].notna().sum()}")
    print(f"Database matches (SV): {len(sv_annotated[sv_annotated['Common_name'] != 'NOT IN CATALOGUE'])}")
    print(f"\nOutput directory: {output_dir}")
    

