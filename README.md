# NanoThal

ONT variant calling for α/β-thalassaemia. Targets HBA1/HBA2 (chr16) and HBB (chr11).

Clair3 (SNV/indel, phased) + Sniffles and CuteSV (SV) → region filter → VEP annotation
→ clinical report against a curated variant database (70 variants, UAE cohort).

## Status: mid-development

Validation is currently simulation-based (Badread, `nanopore2023` error model).
`scripts/simulate/make_haplotype.py` splices published breakpoints out of hg38 and
emits a truth VCF; `simulation/truth_a37.vcf` is the -α3.7 answer key.

Results so far, on simulated data with known ground truth:
- Wild-type null control: CuteSV 38 false positives at min_support 3, 0 at 25.
  Sniffles 0 throughout.
- Het -α3.7 across 930x → 30x: support falls 398 → 10, fraction holds at ~40%.
  Sniffles calls 0/1 correctly at every depth; size exact (3812bp) at every depth.
- Implication: a fixed min_support count cannot generalise across 1,300x amplicon
  and 30x WGS. Threshold should be a fraction of median depth.

## Known issues

- SV matching (`scripts/identify_sv.py`) uses containment, not reciprocal overlap.
  Produces false positives, e.g. "Hb Lepore" (a ~7kb deletion) called from a 92bp deletion.
- `cutesv_min_support: 3` is known-wrong; kept until threshold work lands.
- CuteSV reports `./.` even with `--genotype` and index present (DV counted, DR not).
  Zygosity currently taken from Sniffles.
- VEP REST returns `unknown` for every `upstream_gene_variant` — i.e. exactly where the
  β-thal promoter variants sit (c.-151C>T, c.-138C>A). Indels return api_error.
- Co-inheritance flag fires on any HBA+HBB variant, not just causative ones.
- `annotate_structural` handles HBA only; HBB structural variants are never annotated.
- WGS mode untested: sniffles/cutesv/clair3 hardcode the amplicon reference.
- Target regions hardcoded (chr16:170000-178000, chr11:5225000-5228000).
- Dependencies split across conda and BlueBEAR modules — Clair3's PYTHONPATH shadows
  conda's pysam, so SV callers and Clair3 cannot share a shell (hence env_sv.sh /
  env_clair3.sh). `--use-conda` is the fix; not yet done.

## Not reproducible outside BlueBEAR yet

Paths in `config.yml` are absolute and local.
