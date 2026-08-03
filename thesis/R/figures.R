source("R/setup.R")

plot_reported_spectrum <- function(data) {
  data |>
    mutate(category = reorder(category, count)) |>
    ggplot(aes(category, count)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(labels = comma) +
    labs(
      x = NULL,
      y = "Reported records",
      title = "Reported diagnostic spectrum",
      caption = "Replace synthetic example values with governed cohort aggregates."
    )
}

plot_assay_coverage <- function(data) {
  data |>
    mutate(
      observable = factor(observable, levels = c(FALSE, TRUE)),
      haplotype_id = factor(haplotype_id, levels = rev(unique(haplotype_id)))
    ) |>
    ggplot(aes(product_id, haplotype_id, fill = observable)) +
    geom_tile(linewidth = 0.2) +
    scale_fill_manual(values = c("grey90", "grey25"), labels = c("No", "Yes")) +
    labs(
      x = "Declared assay product",
      y = "Chromosome haplotype",
      fill = "Observable",
      title = "Assay coverage matrix"
    ) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

plot_molecule_status <- function(data) {
  data |>
    mutate(fraction = count / sum(count), status = reorder(status, fraction)) |>
    ggplot(aes(status, fraction)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = NULL,
      y = "Fraction of reads",
      title = "Primer-aware molecule admission"
    )
}

plot_genotype_posterior <- function(data, n = 10) {
  data |>
    slice_min(rank, n = n) |>
    mutate(genotype_pairs = reorder(genotype_pairs, posterior)) |>
    ggplot(aes(genotype_pairs, posterior)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, 1)) +
    labs(
      x = "Candidate genotype class",
      y = "Conditional research posterior",
      title = "Ranked assay-observable genotype classes"
    )
}

plot_product_efficiency <- function(data) {
  data |>
    mutate(product_id = reorder(product_id, relative_efficiency)) |>
    ggplot(aes(product_id, relative_efficiency)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = 2) +
    coord_flip() +
    labs(
      x = "Compiled PCR product",
      y = "Relative fitted efficiency",
      title = "Control-fitted product efficiencies",
      caption = paste(
        "Efficiencies are relative within identifiable product components;",
        "a single product has no independently identified scale."
      )
    )
}

plot_calibration_fit <- function(data) {
  data |>
    filter(eligible, product_id != "") |>
    ggplot(aes(fitted_weight, observed_weight, shape = stratum)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(size = 2.2) +
    facet_wrap(vars(product_id), scales = "free") +
    labs(
      x = "Fitted effective product weight",
      y = "Observed effective product weight",
      shape = "Calibration stratum",
      title = "Observed versus fitted control-product evidence"
    )
}

plot_backend_timings <- function(data) {
  data |>
    mutate(backend = factor(backend, levels = unique(backend))) |>
    ggplot(aes(backend, seconds)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.06), alpha = 0.7) +
    scale_y_continuous(labels = label_number(accuracy = 0.001)) +
    labs(
      x = NULL,
      y = "Wall-clock seconds",
      title = "Reference Python and grouped Cython scoring runtimes",
      caption = "Replace synthetic timings with repeated measurements on frozen inputs."
    )
}
