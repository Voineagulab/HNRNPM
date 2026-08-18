#!/usr/bin/env bash
# rMATS-turbo: HnrnpM vs NEG4 on PSD5-only samples (timepoint sensitivity).
# 2 HnrnpM_PSD5 vs 2 NEG4_PSD5 — early-KD splicing effects.
# Run on rna2 server. Activate the rmats conda env first.
set -euo pipefail

ROOT=/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication
GTF=/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/DATA/Homo_sapiens.GRCh38.109.gtf

OD=$ROOT/TimepointAnalysis/RESULTS/rmats_PSD5_only_HnrnpM
TMP=$ROOT/TimepointAnalysis/RESULTS/rmats_PSD5_only_tmp
LOG=$ROOT/TimepointAnalysis/RESULTS/rmats_PSD5_only_run.log

mkdir -p "$OD" "$TMP"

source /home/rna2/miniconda3/etc/profile.d/conda.sh
conda activate rmats

rmats.py \
    --b1 "$ROOT/TimepointAnalysis/SCRIPTS/b1_HnrnpM_PSD5_only.txt" \
    --b2 "$ROOT/TimepointAnalysis/SCRIPTS/b2_NEG4_PSD5_only.txt" \
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

echo "rMATS PSD5-only finished. Output in $OD"
