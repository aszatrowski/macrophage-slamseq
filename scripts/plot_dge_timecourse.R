# Activate renv environment
source("renv/activate.R")
library(ggplot2)
suppressWarnings(library(dplyr))

summary_stats <- readr::read_csv(snakemake@input$dge_summary_stats, show_col_types = FALSE)
fdr_threshold <- snakemake@params$fdr_threshold
logFC_threshold <- snakemake@params$logFC_threshold
n_timepoints_dge_threshold <- snakemake@params$n_timepoints_dge_threshold
max_genes_to_plot <- snakemake@params$max_genes_to_plot 

# snakemake provides params twice, once named and once unnamed; extract only the named ones and construct a string for plot subtitle
keys <- names(snakemake@params)
keys <- keys[!keys == ""]
values <- unlist(snakemake@params[keys])
params_string <- paste(paste(keys, values, sep = ": "), collapse = ", ")

# function: extract all genes passing fdr_threshold% FDR and |FC| > FC_threshold
extract_de_genes <- function(df, fdr_threshold = 0.05, logFC_threshold = 1) {
  de_genes <- dplyr::filter(df, FDR <= fdr_threshold & abs(logFC) >= logFC_threshold)
  return(de_genes)
}

# use by() to call extract_de_genes() over all values of "timepoint", return a named list of DFs with summary stats
de_genes_list <- by(summary_stats, summary_stats$timepoint, extract_de_genes)

# retrieve genes that are DE in >50% of conditions, up to max_genes_to_plot
de_genes_all <- do.call(rbind, de_genes_list)
de_gene_counts <- sort(table(de_genes_all$Gene), decreasing = TRUE)
de_genes_frequent <- head(
  names(de_gene_counts[de_gene_counts > n_timepoints_dge_threshold]),
  max_genes_to_plot
)

# create a DF with reference FC=0 at timepoint 0 for all genes in de_genes_frequent
reference_rows_df <- data.frame(matrix(
  nrow = length(de_genes_frequent),
  ncol = length(colnames(summary_stats))
))
colnames(reference_rows_df ) <- colnames(summary_stats)
reference_rows_df$Gene <- de_genes_frequent
reference_rows_df$timepoint <- 0
reference_rows_df$logFC <- 0

# subset summary_stats to only include genes in de_genes_frequent
de_genes_timecourse <- dplyr::filter(summary_stats, Gene %in% de_genes_frequent)
# append reference rows to de_genes_timecourse so they'll show up in the plot
de_genes_timecourse <- rbind(de_genes_timecourse, reference_rows_df)

# Manual port of palette "Klein" from package MoMAColors 
palette <- c("#FF4D6FFF", "#579EA4FF", "#DF7713FF", "#F9C000FF", "#86AD34FF", "#5D7298FF", "#81B28DFF", "#7E1A2FFF", "#2D2651FF", "#C8350DFF", "#BD777AFF")

timecourse_dge_plot <- ggplot(de_genes_timecourse, aes(x = timepoint, y = logFC, color = Gene)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = rep(palette, length.out = length(de_genes_frequent))) +
  scale_x_continuous(breaks = unique(de_genes_timecourse$timepoint)) +
  theme_minimal() +
  labs(
    x = paste0("Timepoint ", "(", snakemake@params$time_unit, ")"),
    y = bquote(log[2] ~ "(FC)"),
    title = paste0(
      "Time series for DEGs in >",
      n_timepoints_dge_threshold,
      " timepoints, ",
      snakemake@wildcards$readtype,
      " RNA"
    ),
    subtitle = params_string
  )

ggsave(
  filename = snakemake@output$timecourse_dge_plot,
  plot = timecourse_dge_plot,
  width = 10,
  height = 6,
)