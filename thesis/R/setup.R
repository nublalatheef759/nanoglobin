suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(knitr)
  library(readr)
  library(scales)
  library(tibble)
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

require_input <- function(path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(normalized)) {
    stop("Required analysis input is missing: ", path, call. = FALSE)
  }
  normalized
}

read_required_tsv <- function(path, ...) {
  readr::read_tsv(require_input(path), show_col_types = FALSE, ...)
}

read_required_csv <- function(path, ...) {
  readr::read_csv(require_input(path), show_col_types = FALSE, ...)
}

read_required_json <- function(path, ...) {
  jsonlite::read_json(require_input(path), simplifyVector = TRUE, ...)
}
