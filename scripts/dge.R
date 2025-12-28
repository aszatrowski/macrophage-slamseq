library(edgeR)
read_counts_df <- readr::read_csv(snakemake@input$merged_counts, show_col_types = FALSE) |>
  dplyr::mutate(across(where(is.double), ~ replace(.x, is.na(.x), 0)))

timepoints <- as.factor(gsub("t_(.+)m_(nascent|total).*", "\\1", colnames(read_counts_df)[-1]))

dge <- DGEList(
  counts = read_counts_df[, -1],
  group = timepoints,
  genes = read_counts_df$Gene
)

# Filter lowly expressed genes. additional min.counts are required for NAs/too-low counts in total set
# Require at least min.count reads in at least 2 samples (one replicate pair)
keep <- filterByExpr(dge, min.count = 10, min.total.count = 15)
# Additional: remove genes with zeros in more than 50% of samples
keep2 <- rowSums(dge$counts > 0) >= 7  # At least 7/14 samples non-zero

dge <- dge[(keep & keep2), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix for timepoints: compares each timepoint to baseline (t_0m), with intercepts, no interactions
design <- model.matrix(~ timepoints, data = dge$samples)

fit <- glmQLFit(dge, design)


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
  results <- results |>
    dplyr::left_join(symbol_ensg_mapping, by = c("genes" = "Gene")) |>
    dplyr::rename(Gene = "Symbol", ENSG = "genes") |>
    dplyr::mutate(
      timepoint = timepoint,
      comparison = paste(timepoint, "vs", "0m", sep = "_")) |>
    dplyr::select(Gene, timepoint, comparison, logFC, logCPM, F, PValue, FDR, ENSG)
  return(results)
}

symbol_ensg_mapping <- readr::read_csv(snakemake@input$symbol_ensg_mapping, show_col_types = FALSE)

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
readr::write_csv(dge_results, snakemake@output$dge_summary_stats)