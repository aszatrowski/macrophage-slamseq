library(tidyverse)
tsv <- read_tsv("data/slam_quant/s27_no4sU_test/grandslam.tsv") |>
  rename_all(~ str_replace(., "LB-HT-28s-JL-09_S27", "sample")) |>
  rename_all(~ str_replace(., "LB-HT-28s-HT-09_S9.no4sU", "control")) |>
  rename_all(~ str_replace(., " ", "_"))

map_histogram <- ggplot(tsv, aes(x = sample_MAP)) +
  geom_histogram()

ggsave(
  "outputs/figures/test_gs_results_sample_map_histogram.png",
  plot = map_histogram,
  width = 6,
  height = 4,
  dpi = 300
)
plot_nascent_total_bygene <- function(df, gene_ids) {
  df_genes <- df |>
    filter(Symbol %in% gene_ids) |>
    mutate(
      nascent_count = sample_Readcount * sample_MAP,
      old_count = sample_Readcount * (1 - sample_MAP)
    ) |>
    select(Symbol, old_count, nascent_count) |>
    pivot_longer(
      cols = c(old_count, nascent_count),
      names_to = "Measurement",
      values_to = "Count"
    )

  g <- ggplot(df_genes, aes(x = Symbol, y = Count)) +
    geom_bar(aes(fill = Measurement), stat = "identity", position = "stack") +
    scale_y_log10() +
    labs(
      title = paste("Nascent & old counts by Gene"),
      x = "gene",
      y = "nascent reads, log scale"
    ) +
    theme_bw()
  ggsave(
    paste0("outputs/figures/nascent_counts_", paste(gene_ids, collapse = "_"), ".png"),
    plot = g,
    width = 6,
    height = 4,
  )
}
plot_nascent_total_bygene(tsv, c("GAPDH", "MYC", "FOS", "JUN", "TNF", "IL1", "IL6", "IL8", "NFKB1", "TGFB1"))