library(ggplot2)
fdr_threshold <- snakemake@params$fdr_threshold
logFC_threshold <- snakemake@params$logFC_threshold
comparison <- snakemake@wildcards$comparison

plot_width <- snakemake@params$plot_width
plot_height <- snakemake@params$plot_height

summary_stats <- readr::read_csv(snakemake@input$dge_summary_stats, show_col_types = FALSE) |>
  dplyr::filter(comparison == !!comparison) |>
  dplyr::mutate(
    significance = dplyr::case_when(
      FDR < fdr_threshold & logFC > logFC_threshold ~ "Upregulated",
      FDR < fdr_threshold & logFC < -logFC_threshold ~ "Downregulated",
      .default = "Not Significant"
    )
  )

volcano_plot <- ggplot(summary_stats, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey"
    )
  ) +
  ggrepel::geom_text_repel(
    data = subset(summary_stats, FDR < fdr_threshold & abs(logFC) > logFC_threshold),
    aes(label = Gene),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    x = bquote(log[2] ~ "(FC)"),
    y = bquote(-log[10] ~ "(FDR)"),
    title = paste("Volcano Plot:", comparison),
  ) +
  theme(legend.title = element_blank())

ggsave(
  filename = snakemake@output$volcano_plot,
  plot = volcano_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
)