read_counts_list <- lapply(snakemake@input$read_counts, readr::read_csv, show_col_types = FALSE)
# Remove columns that are not shared across all data frames
shared_columns <- intersect(
  names(read_counts_list[[1]]),
  names(read_counts_list[[2]])
)
read_counts_list <- lapply(read_counts_list, function(x) dplyr::select(x, all_of(shared_columns)))
merged_counts_df <- purrr::reduce(
  read_counts_list,
  function(x, y) dplyr::full_join(
    x, y,
    by = "Gene",
    suffix = paste0("_", snakemake@params$donor)
  )
)
readr::write_csv(merged_counts_df, snakemake@output$merged_counts )