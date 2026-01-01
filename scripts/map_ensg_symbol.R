# Activate renv environment
source("renv/activate.R")
read_counts <- readr::read_tsv(snakemake@input[[1]], show_col_types = FALSE)
symbol_ensg_mapping <- read_counts |>
  dplyr::select("Gene", "Symbol") |>
  dplyr::distinct()
readr::write_csv(symbol_ensg_mapping, snakemake@output$symbol_ensg_mapping)