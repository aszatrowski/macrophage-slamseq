library(edgeR)
read_counts <- readr::read_csv(snakemake@input$merged_counts, show_col_types = FALSE) |>
  dplyr::mutate(across(where(is.double), ~ replace(.x, is.na(.x), 0)))

counts_matrix <- as.matrix(read_counts[, -1])

metadata <- readr::read_csv(snakemake@input$metadata, show_col_types = FALSE)

timepoint <- factor(metadata$timepoint)
donor <- factor(metadata$donor)
replicate <- factor(metadata$replicate)

dge <- DGEList(
  counts = counts_matrix,
  samples = metadata,
  genes = read_counts[, c("Gene")]
)
dge <- calcNormFactors(dge, method = "TMM")

# create design matrix; reference sample is donor 1, timepoint 0. if there are multiple donors, include donor in the design matrix
if (length(unique(donor)) > 1) {
  design <- model.matrix(~ timepoint + replicate + donor, data = metadata)
} else {
  design <- model.matrix(~ timepoint + replicate, data = metadata)
}

dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)

qlf_t1_t0 <- glmQLFTest(fit, coef = 2) 
tt <- topTags(qlf_t1_t0, n = 50)
symbol_ensg_mapping <- readr::read_csv(snakemake@input$symbol_ensg_mapping, show_col_types = FALSE)
summstats <- tt$table |>
  dplyr::left_join(symbol_ensg_mapping, by = c("Gene" = "Gene")) |>
  dplyr::rename(Gene = "Symbol", ENSG = "Gene") |>
  dplyr::select(Gene, ENSG, logFC, logCPM, F, PValue, FDR)
readr::write_csv(summstats, snakemake@output$dge_summary_stats)