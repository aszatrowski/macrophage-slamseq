library(ggplot2)
set.seed(100) # for reproducibility of ggrepel label positions; otherwise unnecessary

fdr_threshold <- snakemake@params$fdr_threshold
logFC_threshold <- snakemake@params$logFC_threshold
comparison <- snakemake@wildcards$comparison

palette <- snakemake@params$palette

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
  geom_point(aes(color = significance), size = 1.5) +
  scale_color_manual(values = c("Upregulated" = palette[[2]], "Downregulated" = palette[[1]], "Not Significant" = palette[[4]])) +
  ggrepel::geom_label_repel(
    data = subset(summary_stats, FDR < fdr_threshold & abs(logFC) > logFC_threshold),
    aes(label = Gene),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.25,
    max.overlaps = 15
  ) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    x = bquote(log[2] ~ "(FC)"),
    y = bquote(-log[10] ~ "(FDR)"),
    title = paste("Volcano Plot:", snakemake@wildcards$readtype, comparison),
  ) +
  theme(legend.title = element_blank())

ggsave(
  filename = snakemake@output$volcano_plot,
  plot = volcano_plot,
  width = plot_width,
  height = plot_height,
)