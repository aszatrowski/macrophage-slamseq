# Activate renv environment
source("renv/activate.R")
suppressWarnings(library(dplyr))
library(ggplot2)

read_summary_stats <- function(path, comparison_value) {
  # regex to retrieve readtype from filename
  readtype <- gsub(".*summary_stats_(nascent|total)\\.csv$", "\\1", path)
  logFC_colname <- paste0("logFC_", readtype)
  summary_stats_readtype <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::filter(comparison == comparison_value) |>
    dplyr::rename(!!logFC_colname := logFC) |>
    select(Gene, !!logFC_colname, ENSG)
  return(summary_stats_readtype)
}

effect_sizes_nascent_total_list <- lapply(
  X = snakemake@input$dge_summary_stats,
  FUN = read_summary_stats,
  comparison = snakemake@wildcards$comparison
)
effect_sizes_nascent_total_df <- Reduce(
  function(x, y) dplyr::inner_join(x, y, by = c("ENSG", "Gene")), # always join by ENSG to avoid non-unique gene names
  effect_sizes_nascent_total_list
) 
palette <- snakemake@params$palette
effect_size_corr_plot <- ggplot(
  effect_sizes_nascent_total_df,
  aes(
    x = logFC_total,
    y = logFC_nascent
  )
) +
  geom_point(alpha = 0.6, color = palette[1]) +
  geom_smooth(method = "lm", color = palette[2], se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, color = palette[4]) +
  coord_fixed() + # equal scaling for x and y axes, so y = x line is at 45 degrees
  labs(
    x = bquote("log"[2] ~ "Fold Change (Total)"),
    y = bquote("log"[2] ~ "Fold Change (Nascent)"),
    title = paste("Effect Size Correlation:", snakemake@wildcards$comparison)
  ) +
  theme_bw()

ggsave(
  filename = snakemake@output$corr_plot,
  plot = effect_size_corr_plot,
  width = 4,
  height = 8
)