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
