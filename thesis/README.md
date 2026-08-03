# Thesis workspace

This is a Quarto book rendered to HTML and Typst PDF from the same source.

```bash
cd thesis
quarto preview                       # interactive HTML review
quarto render --to html              # HTML book
quarto render --to typst             # PDF via bundled Typst
```

The Typst render uses `keep-typ: true`, so the generated `.typ` source remains in
`_book/` for final layout inspection. Quarto bundles the Typst CLI.

## Data policy

No patient-level clinical data or reports are committed. Chapters look first for
approved generated results under `../results/`; if those files are absent, the
book renders against clearly labelled synthetic examples under `data/example/`.
Replace examples only through governed pipeline outputs, not by editing figures.

Synthetic fallback tables are publication scaffolding, not reported study results.
Every final figure must identify its generating command, repository commit,
configuration and source-data governance status.

## Figure inputs

Cohort and assay mechanics:

- `results/thesis/cohort_reported_spectrum.csv`
- `results/assay/coverage.tsv`
- `results/thesis/molecule_status.csv`
- `results/assay/genotypes/<sample>.posteriors.tsv`

Calibration and model execution:

- `results/assay/calibration/product_parameters.tsv`
- `results/assay/calibration/sample_product_fit.tsv`
- `results/assay/calibration/calibration.json`
- `results/assay/genotypes_calibrated/<sample>.posteriors.tsv`
- `results/benchmarks/genotype_backend_timings.csv`

Truth scope:

- `results/truth/<benchmark>.globin_coverage.tsv`
- exact benchmark VCF/BED release identifiers and checksums in the software
  appendix or result manifest.

The R helpers in `R/` centralise themes, input resolution and plot definitions.
They include cohort-spectrum, assay-coverage, molecule-QC, product-efficiency,
calibration-residual, posterior and backend-runtime figures.

## Reproducible calibration and evaluation

The thesis must keep the following sets separate:

1. assay-development and synthetic falsification data;
2. control data used to estimate product efficiencies and dropout;
3. parameter-sensitivity or lot/run diagnostics; and
4. the final held-out blinded evaluation set.

The model configuration and calibration JSON are frozen before evaluating the
held-out set. DRAGEN or vendor Reporter output is labelled comparator evidence
unless an event was independently adjudicated.

The validation acquisition registry and event-specific truth requirements are in
[`../validation/README.md`](../validation/README.md) and
[`../validation/resources.tsv`](../validation/resources.tsv).
