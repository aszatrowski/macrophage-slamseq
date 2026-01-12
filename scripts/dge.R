# Activate renv environment
source("renv/activate.R")
library(org.Hs.eg.db)
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
export_results <- function(qlf, timepoint) {
  results <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
  results <- results |>
    dplyr::mutate(
      # rename "genes" with ENSGs to "ENSG"
      ENSG_ei = genes,
      # retrieve Gene column gene symbols, preserving intronic/exonic suffixes; temporarily remove suffixes to match in mapping
      Gene = {
        # remove intronic/exonic suffix for mapping and retrieve symbol from org.Hs.eg.db
        symbols <- mapIds(
          org.Hs.eg.db,
          keys = sub("_(intronic|exonic)$", "", genes),
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
        # Use ENSEMBL ID when symbol is NA
        symbols <- ifelse(is.na(symbols), 
                          sub("_(intronic|exonic)$", "", genes), 
                          symbols)
        # append intronic/exonic suffix back to symbols
        paste(
          symbols,
          dplyr::case_when(
            stringr::str_detect(ENSG_ei, "_intronic$") ~ "intronic",
            stringr::str_detect(ENSG_ei, "_exonic$") ~ "exonic",
            TRUE ~ ""
          ),
          sep = "_"
        )
      },
      timepoint = timepoint,
      comparison = paste(timepoint, "vs", "0m", sep = "_")) |>
    dplyr::select(Gene, timepoint, comparison, logFC, logCPM, F, PValue, FDR, ENSG_ei)
  return(results)
}

dge_results <- do.call(
  rbind,
  mapply(
    export_results,
    qlf = qlf_list,
    timepoint = names(qlf_list),
    SIMPLIFY = FALSE
  )
)
readr::write_csv(dge_results, snakemake@output$dge_summary_stats)