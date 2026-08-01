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

## Figure inputs

- `results/thesis/cohort_reported_spectrum.csv`
- `results/assay/coverage.tsv`
- `results/thesis/molecule_status.csv`
- `results/assay/genotypes/<sample>.posteriors.tsv`

The R helpers in `R/` centralise themes, input resolution and plot definitions.
