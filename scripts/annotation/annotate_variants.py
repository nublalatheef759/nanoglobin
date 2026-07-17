"""
nano-thal: Clinical annotation of thalassaemia variants
Matches pipeline output against curated variant database
and generates clinical interpretation report.
"""

import pandas as pd
import sys
import os
import re
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

def annotate_structural(patient_summary_df, db):
    """Match structural variants (deletions, duplications) against the database"""
    results = []
    
    sv_lookup = {
        '-a3.7': '-α3.7',
        '-a4.2': '-α4.2',
        'anti-3.7': 'ααα(anti-3.7)',
        'anti-4.2': 'ααα(anti-4.2)',
        '--SEA': '--SEA',
        '--FIL': '--FIL',
        '--MED': '--MED',
        'Hb Lepore': 'Hb Lepore'
    }
    
    for _, row in patient_summary_df.iterrows():
        sample = row['Sample']
        
        # Check HBA
        if pd.notna(row.get('HBA')):
            hba_variants = str(row['HBA']).split('; ')
            hba_zygosities = str(row.get('HBA_zygosity', '')).split('; ')
            
            for i, var in enumerate(hba_variants):
                var = var.strip()
                zyg = hba_zygosities[i].strip() if i < len(hba_zygosities) else 'unknown'
                
                # Look up in database
                db_name = sv_lookup.get(var, var)
                db_match = db[db['HVGS'] == db_name]
                
                if not db_match.empty:
                    match = db_match.iloc[0]
                    results.append({
                        'Sample': sample,
                        'Variant_pipeline': var,
                        'Variant_db': db_name,
                        'Common_name': match.get('Common name', ''),
                        'Gene': 'HBA',
                        'Zygosity': zyg,
                        'Type': 'Structural',
                        'Functionality': match.get('Functionality', ''),
                        'In_IthaGenes': match.get('In Ithagenes? (Yes/No)', ''),
                        'In_gnomAD': match.get('In gnomAD? (Yes/No)', ''),
                        'gnomAD_AF': match.get('gnomAD total AF (%)', ''),
                        'gnomAD_ME_AF': match.get('gnomAD ME AF (%)', ''),
                        'Cohort_frequency': match.get('Cohort frequency (%)', '')
                    })
                else:
                    results.append({
                        'Sample': sample,
                        'Variant_pipeline': var,
                        'Variant_db': db_name,
                        'Common_name': 'NOT IN DATABASE',
                        'Gene': 'HBA',
                        'Zygosity': zyg,
                        'Type': 'Structural',
                        'Functionality': 'Unknown',
                        'In_IthaGenes': 'No',
                        'In_gnomAD': 'No',
                        'gnomAD_AF': '',
                        'gnomAD_ME_AF': '',
                        'Cohort_frequency': ''
                    })
    
    return pd.DataFrame(results)

def check_coinheritance(patient_summary_df):
    """Flag patients with both HBA and HBB variants"""
    flags = []
    
    for _, row in patient_summary_df.iterrows():
        has_hba = pd.notna(row.get('HBA')) and row['HBA'] != ''
        has_hbb = pd.notna(row.get('HBB')) and row['HBB'] != ''
        
        if has_hba and has_hbb:
            flags.append({
                'Sample': row['Sample'],
                'HBA': row.get('HBA', ''),
                'HBA_zygosity': row.get('HBA_zygosity', ''),
                'HBB': row.get('HBB', ''),
                'HBB_zygosity': row.get('HBB_zygosity', ''),
                'Flag': 'CO-INHERITED: Both HBA and HBB variants detected. '
                        'HbA2 may be reduced below diagnostic cutoff. '
                        'Standard HPLC screening may miss carrier status.',
                'Risk_level': 'HIGH' if 'homozygous' in str(row.get('HBA_zygosity', '')).lower() 
                             else 'MODERATE'
            })
    
    return pd.DataFrame(flags)

def classify_mutation_type(hgvs, gene):
    """Classify HBB variants by mutation type"""
    if gene != 'HBB' and gene != 'HBB':
        return 'N/A'
    
    beta0 = ['c.93-22_95del', 'c.92+1G>A', 'c.92+1G>T', 'c.118C>T',
             'c.25_26delAA', 'c.17_18delCT', 'c.315+1G>A', 
             'c.27dupG', 'c.126_129delCTTT', 'c.112delT']
    severe_beta_plus = ['c.92+5G>C', 'c.93-21G>A']
    mild_beta_plus = ['c.-151C>T', 'c.-138C>A', 'c.-50A>C', 'c.-137C>G']
    structural = ['c.20A>T', 'c.364G>C', 'c.79G>A', 'c.364G>A']
    
    if pd.isna(hgvs):
        return 'Unknown'
    
    c_part = hgvs.split(':')[1] if ':' in str(hgvs) else str(hgvs)
    
    if c_part in beta0:
        return 'β⁰ (null)'
    elif c_part in severe_beta_plus:
        return 'Severe β⁺'
    elif c_part in mild_beta_plus:
        return 'Mild β⁺'
    elif c_part in structural:
        return 'Structural'
    else:
        return 'Unclassified'

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
            report.append(f"  In IthaGenes: {row['In_IthaGenes']}")
            report.append(f"  In gnomAD: {row['In_gnomAD']}")
            if pd.notna(row['gnomAD_AF']) and str(row['gnomAD_AF']).strip() not in ('', 'nan'):
                report.append(f"  gnomAD AF (total): {row['gnomAD_AF']}%")
                report.append(f"  gnomAD AF (ME): {row['gnomAD_ME_AF']}%")
            if pd.notna(row['Cohort_frequency']) and str(row['Cohort_frequency']).strip() not in ('', 'nan'):
                report.append(f"  UAE cohort frequency: {row['Cohort_frequency']}%")
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
    sample_flags = coinheritance_flags[coinheritance_flags['Sample'] == sample]
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
    db_path = os.path.join(base_dir, "databases", "variant_lookup.csv")
    output_dir = os.path.join(base_dir, "results", "clinical_reports")
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("nano-thal Clinical Annotation Pipeline")
    print("=" * 40)
    
    # Load data
    print("Loading variant data...")
    variants = pd.read_csv(variants_path)
    patients = pd.read_csv(patient_path)
    db = load_database(db_path)
    
    print(f"  Samples: {len(patients)}")
    print(f"  SNV/Indel variants: {len(variants)}")
    print(f"  Database variants: {len(db)}")
    
    # Annotate
    print("\nAnnotating SNV/indels...")
    snv_annotated = annotate_snv_indels(variants, db)
    
    print("Annotating structural variants...")
    sv_annotated = annotate_structural(patients, db)
    
    print("Checking co-inheritance...")
    coinheritance = check_coinheritance(patients)
    
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
    print(f"Database matches (SV): {len(sv_annotated[sv_annotated['Common_name'] != 'NOT IN DATABASE'])}")
    print(f"\nOutput directory: {output_dir}")
    

