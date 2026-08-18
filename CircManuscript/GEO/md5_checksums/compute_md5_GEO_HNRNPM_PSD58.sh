#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# compute_md5_GEO_HNRNPM_PSD58.sh
#
# Purpose:
#   Compute MD5 checksums for every file in the HNRNPM 7-sample GEO/SRA
#   submission (14 raw FASTQs + 3 processed data matrices), so the values
#   can be pasted into the "MD5 Checksums" tab of the GEO metadata template,
#   and so the checksum set is reproducible by anyone with the data.
#
# Output:
#   Two tab-delimited sections ("file name<TAB>file checksum"), one for RAW
#   FILES and one for PROCESSED DATA FILES, matching the layout of the GEO
#   template's MD5 tab. Printed to stdout and, unless disabled, written to
#   $OUT.
#
# Usage:
#   bash compute_md5_GEO_HNRNPM_PSD58.sh
#   ROOT=/mnt/Scratch/PROJECTS/JuliWang/circRBP_pilot bash compute_md5_GEO_HNRNPM_PSD58.sh
#   OUT=/path/to/md5.tsv bash compute_md5_GEO_HNRNPM_PSD58.sh
#   OUT= bash compute_md5_GEO_HNRNPM_PSD58.sh          # stdout only, no file
#
# Notes:
#   - ROOT defaults to the mounted-share path (JW_Katana). On Katana itself
#     set ROOT=/mnt/Scratch/PROJECTS/JuliWang/circRBP_pilot .
#   - md5sum is read-only; it reads (~60 GB of FASTQ) but writes nothing
#     except the small checksum list.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="${ROOT:-/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot}"

RAW_DIR="$ROOT/DATA/WAN11728_merged"
CIRC_DIR="$ROOT/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE"
PROC_DIR="$ROOT/IV/JW/CircManuscript/GEO/processed_files"
MD5_DIR="$ROOT/IV/JW/CircManuscript/GEO/md5_checksums"

# Default output goes to the dedicated md5_checksums/ folder; override with
# $OUT, or set OUT= (empty) to print to stdout only.
OUT="${OUT-$MD5_DIR/md5_out.tsv}"

# 14 raw FASTQs (paired-end R1/R2 for the 7 manuscript samples)
raw_files=(
  "$RAW_DIR/5_gHnrnpM_PSD5_S20_ME_L000_R1_001.fastq.gz"
  "$RAW_DIR/5_gHnrnpM_PSD5_S20_ME_L000_R2_001.fastq.gz"
  "$RAW_DIR/6_gHnrnpM_PSD5_S21_ME_L000_R1_001.fastq.gz"
  "$RAW_DIR/6_gHnrnpM_PSD5_S21_ME_L000_R2_001.fastq.gz"
  "$RAW_DIR/7_gHnrnpM_PSD8_S22_ME_L000_R1_001.fastq.gz"
  "$RAW_DIR/7_gHnrnpM_PSD8_S22_ME_L000_R2_001.fastq.gz"
  "$RAW_DIR/9_gNEG4_PSD5_S24_ME_L000_R1_001.fastq.gz"
  "$RAW_DIR/9_gNEG4_PSD5_S24_ME_L000_R2_001.fastq.gz"
  "$RAW_DIR/11_gNEG4_PSD8_S3_ME_L000_R1_001.fastq.gz"
  "$RAW_DIR/11_gNEG4_PSD8_S3_ME_L000_R2_001.fastq.gz"
  "$RAW_DIR/12_gNEG4_PSD5_S4_ME_L000_R1_001.fastq.gz"
  "$RAW_DIR/12_gNEG4_PSD5_S4_ME_L000_R2_001.fastq.gz"
  "$RAW_DIR/13_gNEG4_PSD8_S5_ME_L000_R1_001.fastq.gz"
  "$RAW_DIR/13_gNEG4_PSD8_S5_ME_L000_R2_001.fastq.gz"
)

# 3 processed data matrices (all 7 samples)
processed_files=(
  "$PROC_DIR/salmon_gene_count_HNRNPM_PSD58.tsv"
  "$CIRC_DIR/bsj_matrix_CIRIquant_HnrnpM_PSD58.tsv"
  "$CIRC_DIR/fsj_matrix_CIRIquant_HnrnpM_PSD58.tsv"
)

# md5sum prints "<hash>  <path>"; reformat to "<basename>\t<hash>" for the tab.
emit () {
  local f
  for f in "$@"; do
    if [[ ! -f "$f" ]]; then
      echo "ERROR: missing file: $f" >&2
      exit 1
    fi
    md5sum "$f" | awk '{ n=$1; $1=""; sub(/^ +/,""); m=$0; sub(/.*\//,"",m); print m "\t" n }'
  done
}

{
  echo -e "# MD5 checksums for GEO/SRA submission (HNRNPM 7-sample series)"
  echo -e "# ROOT=$ROOT"
  echo -e "\nRAW FILES"
  echo -e "file name\tfile checksum"
  emit "${raw_files[@]}"
  echo -e "\nPROCESSED DATA FILES"
  echo -e "file name\tfile checksum"
  emit "${processed_files[@]}"
} | if [[ -n "${OUT:-}" ]]; then mkdir -p "$(dirname "$OUT")"; tee "$OUT"; else cat; fi

if [[ -n "${OUT:-}" ]]; then
  echo "Wrote: $OUT" >&2
fi
