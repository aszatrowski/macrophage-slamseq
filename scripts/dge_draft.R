library(edgeR)
nascent_counts_df <- readr::read_csv("outputs/readcounts/merged_counts_nascent.csv") |>
  dplyr::select("Gene", starts_with("t_0m"), starts_with("t_120m")) |>
  dplyr::mutate(across(where(is.double), ~ replace(.x, is.na(.x), 0)))

dge <- DGEList(
	counts = nascent_counts_df[, -1],
  group = c(rep("t_0m", 2), rep("t_120m", 2)),
	genes = nascent_counts_df$Gene
)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ group, data = dge$samples)

fit <- glmQLFit(dge, design)

# Test for differential expression (t_120m vs t_0m)
qlf <- glmQLFTest(fit, coef = 2)
# coef=2 tests the t_120m coefficient (the time effect)

# Extract results
results <- topTags(qlf, n = Inf, sort.by = "Pvalue")
head(results)