# laod data
merged_read_counts <- readr::read_csv(snakemake@input$merged_counts)
# Extract sample names (exclude Gene and Symbol columns)
samples <- grep("^t_\\d+m_(nascent|total)_readcount", colnames(merged_read_counts), value = TRUE)

# Create metadata
metadata <- data.frame(
  sample = samples,
  timepoint = factor(sub("t_(\\d+)m_(nascent|total)_readcount.*", "\\1", samples)),
  donor = factor(sub(".*_(donor\\d+)_rep\\d+", "\\1", samples)),
  replicate = factor(sub(".*_rep(\\d+)", "\\1", samples))
)

readr::write_csv(metadata, snakemake@output$metadata)