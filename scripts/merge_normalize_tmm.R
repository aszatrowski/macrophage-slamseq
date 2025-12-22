library(edgeR)

read_counts_list <- lapply(snakemake@input$read_counts, readr::read_csv)
merged_counts_df <- purrr::reduce(
  read_counts_list,
  function(x, y) dplyr::full_join(
    x, y,
    by = "Gene",
    suffix = paste0("_", snakemake@params$donor)
  )
)
print(colnames(merged_counts_df))


# dge <- DGEList(counts = read_counts[, -1], group = read_counts$group)
readr::write_csv(merged_counts_df, snakemake@output$merged_normalized_counts)