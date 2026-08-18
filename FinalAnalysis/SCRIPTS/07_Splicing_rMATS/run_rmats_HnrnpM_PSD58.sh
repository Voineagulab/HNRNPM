#!/usr/bin/env bash
# rMATS-turbo: HnrnpM vs NEG4 on PSD5/PSD8-matched samples (publication version).
# Run on rna2 server. Activate the rmats conda env first.
set -euo pipefail

ROOT=/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication
GTF=/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/DATA/Homo_sapiens.GRCh38.109.gtf

OD=$ROOT/RESULTS/07_Splicing_rMATS/rmats_raw
TMP=$ROOT/RESULTS/07_Splicing_rMATS/rmats_tmp
LOG=$ROOT/RESULTS/07_Splicing_rMATS/rmats_run.log

mkdir -p "$OD" "$TMP"

source /home/rna2/miniconda3/etc/profile.d/conda.sh
conda activate rmats

rmats.py \
    --b1 "$ROOT/SCRIPTS/07_Splicing_rMATS/b1_HnrnpM_PSD58.txt" \
    --b2 "$ROOT/SCRIPTS/07_Splicing_rMATS/b2_NEG4_PSD58.txt" \
    --gtf "$GTF" \
    --od  "$OD" \
    --tmp "$TMP" \
    -t paired \
    --libType fr-firststrand \
    --readLength 100 \
    --variable-read-length \
    --allow-clipping \
    --nthread 8 \
    --cstat 0.0001 \
    2>&1 | tee "$LOG"

echo "rMATS finished. Output in $OD"
