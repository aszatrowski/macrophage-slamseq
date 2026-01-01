# Activate renv environment
source("renv/activate.R")

library(dplyr)
library(VennDiagram)

fdr_threshold <- snakemake@params$fdr_threshold
logFC_threshold <- snakemake@params$logFC_threshold

read_summary_stats_get_degs <- function(path, comparison) {
  deg_set <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::filter(comparison == !!comparison) |>
    dplyr::filter(FDR < fdr_threshold & abs(logFC) > logFC_threshold) |>
    select(ENSG) |>
    unlist()
  return(deg_set)
}

deg_set_list <- lapply(
  X = snakemake@input$dge_summary_stats,
  FUN = read_summary_stats_get_degs,
  comparison = snakemake@wildcards$comparison
)

# Check if any DEG sets are empty; creating a dummy file is necessary since venn.diagram() fails on empty sets
check_no_degs <- sapply(deg_set_list, length) == 0
if (any(check_no_degs)) {
  warning("One or more DEG sets are empty. Venn diagram will not be generated.")
  file.create(snakemake@output$venn_diagram)
  quit(status = 0)
}

VennDiagram::venn.diagram(
  disable.logging = TRUE,
  x = deg_set_list,
  category.names = c("nascent", "total"),
  filename = snakemake@output$venn_diagram,
  output=TRUE,

  fill = snakemake@params$palette[1:length(deg_set_list)],
  alpha = 0.5,
  lwd = 0,
  lty = "blank",

  cex = 10,
  fontface = "bold",
  fontfamily = "sans",

  cat.cex = 10,
  cat.fontface = "bold",
  cat.fontfamily = "sans",

  margin = 0.1,

  # VennDiagram takes sizes in pixels and can't guess file type from filename
  imagetype="png",
  height = 8000,
  width = 8000,
  units = "px",
  resolution = 300
)