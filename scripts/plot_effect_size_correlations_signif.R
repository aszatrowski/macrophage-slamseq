# Activate renv environment
source("renv/activate.R")
suppressWarnings(library(dplyr))
library(ggplot2)

intron_exon <- snakemake@wildcards$intron_exon

read_summary_stats <- function(path, comparison_value) {
  # regex to retrieve readtype from filename
  readtype <- gsub(".*summary_stats_(nascent|total)\\.csv$", "\\1", path)
  logFC_colname <- paste0("logFC_", readtype)
  FDR_colname <- paste0("FDR_", readtype)
  summary_stats_readtype <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::filter(stringr::str_detect(Gene, paste0("_", intron_exon))) |>
    dplyr::filter(comparison == comparison_value) |>
    dplyr::rename(!!logFC_colname := logFC) |>
    dplyr::rename(!!FDR_colname := FDR) |>
    select(Gene, !!logFC_colname,!!FDR_colname, ENSG_ei)
  return(summary_stats_readtype)
}

effect_sizes_nascent_total_list <- lapply(
  X = snakemake@input$dge_summary_stats,
  FUN = read_summary_stats,
  comparison = snakemake@wildcards$comparison
)

significance_deg_set <- snakemake@wildcards$deg_set

logFC_colname <- paste0("logFC_", significance_deg_set)
FDR_colname <- paste0("FDR_", significance_deg_set)

effect_sizes_nascent_total_df <- Reduce(
  function(x, y) dplyr::inner_join(x, y, by = c("ENSG_ei", "Gene")), # always join by ENSG_ei to avoid non-unique gene names
  effect_sizes_nascent_total_list
) |> 
  # filter to only significant DEGs in wildcard-specified DEG set
  dplyr::filter(
    !!as.symbol(FDR_colname) <= as.numeric(snakemake@params$fdr_threshold) &
    abs(!!as.symbol(logFC_colname)) > as.numeric(snakemake@params$logFC_threshold)
  )


# safely handle case where there are no significant DEGs—will return NULL instead of an lm object
effect_size_lm <- tryCatch(
  {
    if (nrow(effect_sizes_nascent_total_df) == 0) {
      NULL
    } else {
      lm(logFC_nascent ~ logFC_total, data = effect_sizes_nascent_total_df)
    }
  },
  error = function(e) NULL
)

palette <- snakemake@params$palette
effect_size_corr_plot <- ggplot(
    effect_sizes_nascent_total_df,
    aes(x = logFC_total, y = logFC_nascent)
  ) +
  geom_point(alpha = 0.6, color = palette[1]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, color = palette[4]) +
  labs(
    x = bquote("log"[2] ~ "Fold Change (Total)"),
    y = bquote("log"[2] ~ "Fold Change (Nascent)"),
    title = paste("FC Correlation:", snakemake@wildcards$comparison, intron_exon, "reads; subset to DEGs in", significance_deg_set, "reads"),
    subtitle = bquote(
      R^2 ~ "=" ~ .(if (!is.null(effect_size_lm)) round(summary(effect_size_lm)$r.squared, 3) else "NA")
    )
  ) +
  theme_bw()

# Add the fitted regression line only if lm succeeded and coefficients are finite
if (!is.null(effect_size_lm) && length(coef(effect_size_lm)) >= 2 && all(is.finite(coef(effect_size_lm)[1:2]))) {
  effect_size_corr_plot <- effect_size_corr_plot +
    geom_abline(
      intercept = coef(effect_size_lm)[1],
      slope = coef(effect_size_lm)[2],
      color = palette[2],
      linewidth = 1
    )
}

ggsave(
  filename = snakemake@output$corr_plot,
  plot = effect_size_corr_plot,
  width = 8,
  height = 8
)