# Activate renv environment
source("renv/activate.R")
library(dplyr, quietly = TRUE)
# at runtime, snakemake adds an environment object with all input and output paths as strings
# snakemake@input$read_table is a string containing the path to the read table file
timepoint_counts <- readr::read_tsv(snakemake@input$read_table, show_col_types = FALSE) |>
  dplyr::filter(!stringr::str_detect(Gene, "_intronic")) |>
  # replace NA readcounts with 0. I'm not sure the biological interpretation of NA here, but NA representing undetected transcript seems reasonable, and low counts will get filtered out by edgeR anyway.
  dplyr::mutate(across(where(is.double), ~ replace(.x, is.na(.x), 0))) |>
  dplyr::select(
    "Gene",
    "Symbol",
    contains("Readcount"),
    # "MAP" is the Bayesian maxium a posteriori estimate of the readcount
    # keep this and discard the rest of the distribution specification
    contains("MAP"),
    -contains("no4sU") # drop no4sU readcount columns
  ) |>
  # compatibility: replace all spaces with underscores
  dplyr::rename_with(~ stringr::str_replace_all(., " ", "_")) |>
  # compatibility: prepend "t_" to all columns that start with a number
  dplyr::rename_with(~ paste0("t_", .x), .cols = matches("^[0-9]"))

# get the timepoint prefixes from the columns; using MAP instead of Readcount correctly excludes the control.
# "$" is used to match the end of the string
timepoint_prefixes <- stringr::str_subset(names(timepoint_counts), "_MAP$") |>
  stringr::str_remove("_MAP$")
colnames(timepoint_counts) <- stringr::str_replace_all(colnames(timepoint_counts), "Readcount", "total_readcount")
# For each timepoint, calculate nascent readcount as total_readcount * MAP
for (prefix in timepoint_prefixes) {
  timepoint_counts <- timepoint_counts |>
    # the somewhat convoluted syntax is needed to create a new column with a dynamic name
    # see documentation on := operator in dplyr::mutate
    dplyr::mutate("{prefix}_nascent_readcount" := .data[[paste0(prefix, "_total_readcount")]] * 
                                        .data[[paste0(prefix, "_MAP")]])
}
total_counts_table <- timepoint_counts |>
  dplyr::select(
    "Gene",
    contains("_total_readcount")
  )
nascent_counts_table <- timepoint_counts |>
  dplyr::select(
    "Gene",
    contains("_nascent_readcount")
  )
readr::write_csv(
  nascent_counts_table,
  snakemake@output$nascent_counts
)
readr::write_csv(
  total_counts_table,
  snakemake@output$total_counts
)