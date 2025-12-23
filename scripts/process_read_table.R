library(dplyr, quietly = TRUE)
# at runtime, snakemake adds an environment object with all input and output paths as strings
# snakemake@input$read_table is a string containing the path to the read table file
timepoint_counts <- readr::read_tsv(snakemake@input$read_table, show_col_types = FALSE) |>
  filter(!stringr::str_detect(Gene, "_intronic")) |>
  dplyr::select(
    "Gene",
    "Symbol",
    contains("Readcount"),
    # "MAP" is the Bayesian maxium a posteriori estimate of the readcount
    # keep this and discard the rest of the distribution specification
    contains("MAP")
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
# Create nascent count columns
for (prefix in timepoint_prefixes) {
  timepoint_counts <- timepoint_counts |>
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