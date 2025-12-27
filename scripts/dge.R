library(edgeR)
nascent_counts_df <- readr::read_csv(snakemake@input$merged_counts, show_col_types = FALSE) |>
  dplyr::mutate(across(where(is.double), ~ replace(.x, is.na(.x), 0)))

timepoints <- as.factor(gsub("t_(.+)m_nascent.*", "\\1", colnames(nascent_counts_df)[-1]))

dge <- DGEList(
  counts = nascent_counts_df[, -1],
  group = timepoints,
  genes = nascent_counts_df$Gene
)


# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Calculate normalization factors for library size and composition bias
dge <- calcNormFactors(dge)

# Design matrix for timepoints: compares each timepoint to baseline (t_0m), with intercepts, no interactions
design <- model.matrix(~ timepoints, data = dge$samples)

fit <- glmQLFit(dge, design)

symbol_ensg_mapping <- readr::read_csv(snakemake@input$symbol_ensg_mapping, show_col_types = FALSE)
# Test for differential expression (t_120m vs t_0m)
qlf_list <- lapply(
  seq(2, length(colnames(design))), # all coefficients except intercept
  function(coef) {
    glmQLFTest(fit, coef = coef)
  }
)
# remove "timepoints" prefix from names
names(qlf_list) <- stringr::str_replace(colnames(design)[-1], "timepoints", "")
# order list by timepoint numeric value
qlf_list <- qlf_list[order(as.numeric(names(qlf_list)))]

# Extract results
export_results <- function(qlf, timepoint, mapping) {
  results <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
  print(head(results))
  results <- results |>
    dplyr::left_join(symbol_ensg_mapping, by = c("genes" = "Gene")) |>
    dplyr::rename(Gene = "Symbol", ENSG = "genes") |>
    dplyr::mutate(
      timepoint = timepoint,
      comparison = paste(timepoint, "vs", "0m", sep = "_")) |>
    dplyr::select(Gene, timepoint, comparison, logFC, logCPM, F, PValue, FDR, ENSG)
  return(results)
}
dge_results <- do.call(
  rbind,
  mapply(
    export_results,
    qlf = qlf_list,
    timepoint = names(qlf_list),
    MoreArgs = list(mapping = symbol_ensg_mapping),
    SIMPLIFY = FALSE
  )
)
head(dge_results)
readr::write_csv(dge_results, snakemake@output$dge_summary_stats)