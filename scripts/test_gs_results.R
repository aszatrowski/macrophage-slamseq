library(tidyverse)
tsv <- read_tsv("data/slam_quant/s27_no4sU_test/grandslam.tsv") |>
  rename_all(~ str_replace(., "LB-HT-28s-JL-09_S27", "sample")) |>
  rename_all(~ str_replace(., "LB-HT-28s-HT-09_S9.no4sU", "control")) |>
  rename_all(~ str_replace(., " ", "_"))
g <- ggplot(tsv, aes(x = sample_MAP)) +
  geom_histogram()

ggsave(
  "outputs/figures/test_gs_results_sample_map_histogram.png",
  plot = g,
  width = 6,
  height = 4,
  dpi = 300
)
