# Activate renv environment
source("renv/activate.R")
suppressWarnings(library(dplyr))
library(ggplot2)

intron_exon <- snakemake@wildcards$intron_exon

read_summary_stats <- function(path, comparison_value) {
  # regex to retrieve readtype from filename
  readtype <- gsub(".*summary_stats_(nascent|total)\\.csv$", "\\1", path)
  logFC_colname <- paste0("logFC_", readtype)
  summary_stats_readtype <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::filter(stringr::str_detect(Gene, paste0("_", intron_exon))) |>
    dplyr::filter(comparison == comparison_value) |>
    dplyr::rename(!!logFC_colname := logFC) |>
    select(Gene, !!logFC_colname, ENSG_ei)
  return(summary_stats_readtype)
}

effect_sizes_nascent_total_list <- lapply(
  X = snakemake@input$dge_summary_stats,
  FUN = read_summary_stats,
  comparison = snakemake@wildcards$comparison
)
effect_sizes_nascent_total_df <- Reduce(
  function(x, y) dplyr::inner_join(x, y, by = c("ENSG_ei", "Gene")), # always join by ENSG to avoid non-unique gene names
  effect_sizes_nascent_total_list
) 

effect_size_lm <- lm(
  logFC_nascent ~ logFC_total,
  data = effect_sizes_nascent_total_df
)

palette <- snakemake@params$palette
effect_size_corr_plot <- ggplot(
    effect_sizes_nascent_total_df,
    aes(x = logFC_total, y = logFC_nascent)
  ) +
  geom_point(alpha = 0.6, color = palette[1]) +
  geom_abline(
    intercept = effect_size_lm$coefficients[1],
    slope = effect_size_lm$coefficients[2],
    color = palette[2],
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, color = palette[4]) +
  labs(
    x = bquote("log"[2] ~ "Fold Change (Total)"),
    y = bquote("log"[2] ~ "Fold Change (Nascent)"),
    title = paste("FC Correlation:", snakemake@wildcards$comparison, intron_exon, "reads"),
    subtitle = bquote("R"^2 ~ "=" ~ .(round(summary(effect_size_lm)$r.squared, 3)))
  ) +
  theme_bw()

ggsave(
  filename = snakemake@output$corr_plot,
  plot = effect_size_corr_plot,
  width = 5,
  height = 8
)