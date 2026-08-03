suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(readr)
  library(scales)
  library(tidyr)
})

options(dplyr.summarise.inform = FALSE, scipen = 999)

theme_nanoglobin <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

theme_set(theme_nanoglobin())

resolve_input <- function(primary, fallback) {
  primary <- normalizePath(primary, mustWork = FALSE)
  fallback <- normalizePath(fallback, mustWork = TRUE)
  if (file.exists(primary)) primary else fallback
}

read_tsv_or_example <- function(primary, fallback, ...) {
  readr::read_tsv(resolve_input(primary, fallback), show_col_types = FALSE, ...)
}

read_csv_or_example <- function(primary, fallback, ...) {
  readr::read_csv(resolve_input(primary, fallback), show_col_types = FALSE, ...)
}

read_json_or_example <- function(primary, fallback, ...) {
  jsonlite::read_json(resolve_input(primary, fallback), simplifyVector = TRUE, ...)
}
