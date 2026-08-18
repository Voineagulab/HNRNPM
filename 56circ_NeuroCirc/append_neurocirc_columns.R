#!/usr/bin/env Rscript
# append_neurocirc_columns.R
# -----------------------------------------------------------------------------
# Append the two NeuroCirc brain-expression flags "Brain1.CTX.CB" and
# "Brain2.DLPFC" to the existing Table S5 in 56_circ.xlsx and save the result
# as 56_circ_neuroscirc.xlsx.
#
# The matching logic is copied from
#   IV/JW/RT_qPCR_validation/build_rt_qpcr_table.R  (step 10, "NeuroCirc flags"):
#     - NeuroCirc IDs are 1-based start, hg38:  chr1_805799_810170_-
#     - Our circRNA IDs are BED-style 0-based start, hg38:  1:805798-810170:-
#       -> shift the start by +1 to build the NeuroCirc key.
#     - Flags are coerced to lowercase "true"/"false" strings (blank if absent).
#
# Input  (read-only):
#   56circ_NeuroCirc/56_circ.xlsx           (existing table; Sheet "S5_BSJxFSJ_circRNA_class")
#   IV/JW/neurocirc.csv                      (NeuroCirc source table)
# Output:
#   56circ_NeuroCirc/56_circ_neuroscirc.xlsx (existing table + 2 appended columns)
#
# Only the output workbook is written; no existing file is modified.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(openxlsx))

HERE  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/56circ_NeuroCirc"
IN_XLSX  <- file.path(HERE, "56_circ.xlsx")
OUT_XLSX <- file.path(HERE, "56_circ_neuroscirc.xlsx")
NEUR     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/neurocirc.csv"

# ----------------------------------------------------------------------------
# 1) Read the existing sheet raw (no header interpretation) so we can locate
#    the header row and the circRNA-ID column exactly as laid out.
#    NB: skipEmptyRows/Cols MUST be FALSE so the reported row/col indices are
#    the true sheet coordinates that loadWorkbook()/writeData() operate on
#    (the default TRUE collapses blank rows and would misalign the writes).
# ----------------------------------------------------------------------------
raw <- read.xlsx(IN_XLSX, sheet = 1, colNames = FALSE,
                 skipEmptyRows = FALSE, skipEmptyCols = FALSE)
cat(sprintf("[1] read %s : %d rows x %d cols\n", basename(IN_XLSX), nrow(raw), ncol(raw)))

# The header row is the one whose first cell is exactly "circRNA ID".
hdr_row <- which(trimws(as.character(raw[[1]])) == "circRNA ID")
stopifnot(length(hdr_row) == 1L)
last_row  <- max(which(!is.na(raw[[1]]) & trimws(as.character(raw[[1]])) != ""))
data_rows <- (hdr_row + 1L):last_row
n_existing_cols <- ncol(raw)
cat(sprintf("[1] header at row %d; %d data rows (rows %d:%d); %d existing cols\n",
            hdr_row, length(data_rows), min(data_rows), max(data_rows), n_existing_cols))

circ_ids <- trimws(as.character(raw[data_rows, 1]))

# ----------------------------------------------------------------------------
# 2) NeuroCirc key builder (identical to build_rt_qpcr_table.R::to_neuro_key)
#    BED-style 0-based "chr:start-end:strand" -> "chrN_(start+1)_end_strand"
# ----------------------------------------------------------------------------
to_neuro_key <- function(circ) {
  if (is.na(circ) || circ == "") return(NA_character_)
  parts <- strsplit(circ, ":", fixed = TRUE)[[1L]]
  stopifnot(length(parts) == 3L)
  chrom  <- parts[1L]
  coords <- strsplit(parts[2L], "-", fixed = TRUE)[[1L]]
  start  <- as.integer(coords[1L])
  end    <- as.integer(coords[2L])
  strand <- parts[3L]
  cc     <- if (startsWith(chrom, "chr")) chrom else paste0("chr", chrom)
  sprintf("%s_%d_%d_%s", cc, start + 1L, end, strand)
}
neuro_key <- vapply(circ_ids, to_neuro_key, character(1))

# ----------------------------------------------------------------------------
# 3) Look up the two flags from neurocirc.csv, coercing to "true"/"false"/"".
# ----------------------------------------------------------------------------
neuro <- read.csv(NEUR, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("circ_id_hg38", "Brain1.CTX.CB", "Brain2.DLPFC") %in% names(neuro)))

idx <- match(neuro_key, neuro$circ_id_hg38)

coerce_flag <- function(v) {
  ifelse(is.na(v), "",
         ifelse(tolower(as.character(v)) %in% c("true", "false"),
                tolower(as.character(v)), ""))
}
brain1 <- coerce_flag(neuro$Brain1.CTX.CB[idx])
brain2 <- coerce_flag(neuro$Brain2.DLPFC[idx])

cat(sprintf("[3] NeuroCirc matched (key found): %d/%d\n",
            sum(!is.na(idx)), length(idx)))
cat(sprintf("[3] Brain1.CTX.CB non-blank: %d ; Brain2.DLPFC non-blank: %d\n",
            sum(brain1 != ""), sum(brain2 != "")))

# ----------------------------------------------------------------------------
# 4) Load the workbook (preserves title, formatting, merged cells) and append
#    the two columns: header cells at hdr_row, values at the data rows.
# ----------------------------------------------------------------------------
wb    <- loadWorkbook(IN_XLSX)
sheet <- names(wb)[1]
c1    <- n_existing_cols + 1L   # Brain1.CTX.CB
c2    <- n_existing_cols + 2L   # Brain2.DLPFC

writeData(wb, sheet, x = "Brain1.CTX.CB", startCol = c1, startRow = hdr_row, colNames = FALSE)
writeData(wb, sheet, x = "Brain2.DLPFC", startCol = c2, startRow = hdr_row, colNames = FALSE)
writeData(wb, sheet, x = brain1, startCol = c1, startRow = min(data_rows), colNames = FALSE)
writeData(wb, sheet, x = brain2, startCol = c2, startRow = min(data_rows), colNames = FALSE)

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("\nWrote: %s  (appended cols %d 'Brain1.CTX.CB' and %d 'Brain2.DLPFC')\n",
            OUT_XLSX, c1, c2))
