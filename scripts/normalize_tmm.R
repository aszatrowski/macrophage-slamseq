library(readr)
library(edgeR)

read_counts <- read_csv(snakemake@input$read_counts)

dge <- DGEList(counts = read_counts[, -1], group = read_counts$group)