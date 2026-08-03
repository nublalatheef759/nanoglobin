source("R/setup.R")

plot_small_variant_performance <- function(data) {
  data |>
    transmute(
      label = paste(region, variant_type, sep = " / "),
      Recall = recall,
      Precision = precision
    ) |>
    pivot_longer(
      cols = c(Recall, Precision),
      names_to = "metric",
      values_to = "value"
    ) |>
    mutate(label = factor(label, levels = rev(unique(label)))) |>
    ggplot(aes(value, label, shape = metric)) +
    geom_point(size = 2.4) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.02)) +
    labs(
      x = NULL,
      y = NULL,
      shape = NULL,
      title = "Truth-scoped small-variant performance"
    )
}
