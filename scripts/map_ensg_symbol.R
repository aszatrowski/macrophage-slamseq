# Activate renv environment
source("renv/activate.R")
read_counts <- readr::read_tsv(snakemake@input[[1]], show_col_types = FALSE)
symbol_ensg_mapping <- read_counts |>
  dplyr::mutate(
    # append "_intronic" or "_exonic" to Gene names for separate identification later
    Symbol = dplyr::case_when(
      # if intronic, leave as is
      stringr::str_detect(Symbol, "_intronic") ~ Symbol,
      # otherwise, append "_exonic"
      TRUE ~ paste0(Symbol, "_exonic")
    ),
    # Gene is ENSG ID
    Gene = dplyr::case_when(
      # if intronic, leave as is
      stringr::str_detect(Symbol, "_intronic") ~ Gene,
      # otherwise, append "_exonic"
      TRUE ~ paste0(Gene, "_exonic")
    )
  ) |>
  dplyr::select("Gene", "Symbol") |>
  dplyr::distinct()
readr::write_csv(symbol_ensg_mapping, snakemake@output$symbol_ensg_mapping)