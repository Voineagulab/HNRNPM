#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# salmon_gene_count_HNRNPM_PSD58.R
#
# Purpose:
#   Generate the GEO/SRA processed data file for the HNRNPM knockdown
#   RNA-seq series by trimming the full 24-sample salmon gene-count matrix
#   down to the seven manuscript samples (3 gHnrnpM + 4 gNEG4 controls,
#   PSD5/PSD8), keeping the gene_id and gene_name annotation columns.
#
#   The source matrix contains all 24 pilot samples and is left UNTOUCHED.
#   A new, separate 7-sample copy is written solely for GEO deposition; no
#   original file and no analysis script is modified.
#
# Input:
#   /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_rnaseq_230423/
#       results/star_salmon/salmon.merged.gene_counts.tsv   (24 samples)
#
# Output:
#   /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/GEO/
#       processed_files/salmon_gene_count_HNRNPM_PSD58.tsv   (7 samples)
#
# Notes:
#   - Column names carry salmon's "X" prefix on the numeric sample IDs.
#   - Katana root here is /mnt/Scratch/PROJECTS/JW_Katana/... (the mounted
#     share). On the Katana host itself the equivalent root is
#     /mnt/Scratch/PROJECTS/JuliWang/... - adjust in_path if running there.
# ---------------------------------------------------------------------------

## ---- Paths -----------------------------------------------------------------
in_path  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_rnaseq_230423/results/star_salmon/salmon.merged.gene_counts.tsv"
out_dir  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/GEO/processed_files"
out_path <- file.path(out_dir, "salmon_gene_count_HNRNPM_PSD58.tsv")

## ---- The seven manuscript samples (order preserved in output) --------------
sample_cols <- c(
  "X5_gHnrnpM_PSD5_S20",   # HNRNPM KD, PSD5
  "X6_gHnrnpM_PSD5_S21",   # HNRNPM KD, PSD5
  "X7_gHnrnpM_PSD8_S22",   # HNRNPM KD, PSD8
  "X9_gNEG4_PSD5_S24",     # Control (gNEG4), PSD5
  "X11_gNEG4_PSD8_S3",     # Control (gNEG4), PSD8
  "X12_gNEG4_PSD5_S4",     # Control (gNEG4), PSD5
  "X13_gNEG4_PSD8_S5"      # Control (gNEG4), PSD8
)
annot_cols <- c("gene_id", "gene_name")

## ---- Read source matrix ----------------------------------------------------
if (!file.exists(in_path)) {
  stop("Source matrix not found: ", in_path)
}
counts <- read.delim(
  in_path,
  header           = TRUE,
  sep              = "\t",
  check.names      = FALSE,   # keep original column names (incl. "X" prefix) verbatim
  stringsAsFactors = FALSE
)

## ---- Validate expected columns ---------------------------------------------
need    <- c(annot_cols, sample_cols)
missing <- setdiff(need, colnames(counts))
if (length(missing) > 0) {
  stop("Expected column(s) not found in source matrix:\n  ",
       paste(missing, collapse = "\n  "),
       "\nAvailable columns:\n  ",
       paste(colnames(counts), collapse = "\n  "))
}

## ---- Subset ----------------------------------------------------------------
out <- counts[, need, drop = FALSE]

message("Source matrix : ", in_path)
message("  dimensions  : ", nrow(counts), " genes x ", ncol(counts), " columns")
message("Trimmed output: ", nrow(out), " genes x ", ncol(out),
        " columns (", length(sample_cols), " samples + ",
        length(annot_cols), " annotation)")

## ---- Write output ----------------------------------------------------------
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
write.table(
  out,
  file      = out_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
message("Wrote: ", out_path)
