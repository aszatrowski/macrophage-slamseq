library(fastqcr)
library(dplyr)
library(readr)

print("starting execute")
fastqc_report <- qc_aggregate(qc.dir = "outputs/fastqc_reports/") |>
    write_csv(snakemake@output[["summary"]])

print("full report written")
fastqc_fails_warnings <- fastqc_report |>
    filter(status != "PASS")

n_warnings <- nrow(filter(fastqc_fails_warnings, status == "WARN"))
n_failures <- nrow(filter(fastqc_fails_warnings, status == "FAIL"))

print(paste(
    n_warnings, "warnings detected."
))
print(paste(
    n_failures, "failures detected."
))

print("writing fails and warnings CSV")
write_csv(fastqc_fails_warnings, snakemake@output[["fails_warnings"]])
print("done")