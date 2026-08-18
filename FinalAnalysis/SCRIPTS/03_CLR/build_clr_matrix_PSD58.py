"""
Build the CLR (circular-to-linear ratio) matrix for the 5,820 CIRIquant-filtered
circRNAs x 7 PSD58 samples.

CLR per (circRNA, sample) is computed directly from the CIRIquant BSJ and FSJ
matrices built in FinalAnalysis step 01:
    CLR = 2*BSJ / (2*BSJ + FSJ)        if (BSJ + FSJ) > 0
        = 0                            otherwise  (no detection in this sample)

This formulation matches CIRIquant's internal `junc_ratio` definition, but
re-derives it directly from the per-sample BSJ + FSJ counts so the CLR matrix
is internally consistent with the count matrices used for BSJ and FSJ DE in
step 01. Universe = 5,820 (same as BSJ DE input).
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS")
BSJ  = ROOT / "01_circRNA_DE" / "bsj_matrix_CIRIquant_HnrnpM_PSD58.tsv"
FSJ  = ROOT / "01_circRNA_DE" / "fsj_matrix_CIRIquant_HnrnpM_PSD58.tsv"
OUT  = ROOT / "03_CLR" / "clr_matrix_HnrnpM_PSD58.tsv"
OUT.parent.mkdir(parents=True, exist_ok=True)

META_COLS = ["circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id"]

bsj = pd.read_csv(BSJ, sep="\t")
fsj = pd.read_csv(FSJ, sep="\t")
samples = [c for c in bsj.columns if c not in META_COLS]
print(f"Samples: {samples}")

bsj_mat = bsj[samples].values.astype(float)
fsj_mat = fsj[samples].values.astype(float)
denom = 2 * bsj_mat + fsj_mat
clr_mat = np.zeros_like(bsj_mat)
mask = denom > 0
clr_mat[mask] = (2 * bsj_mat[mask]) / denom[mask]

out = bsj[META_COLS].copy()
clr_df = pd.DataFrame(clr_mat, columns=samples, index=bsj.index)
out = pd.concat([out, clr_df], axis=1)
out.to_csv(OUT, sep="\t", index=False)

# Summary
print(f"\nCLR matrix shape: {clr_df.shape}")
print(f"  Mean CLR per sample:")
for s in samples:
    grp = "HnrnpM" if "gHnrnpM" in s else "NEG4"
    print(f"    {s:<25} {grp}  mean_CLR = {clr_df[s].mean():.4f}")
group_HM  = [s for s in samples if "gHnrnpM" in s]
group_NEG = [s for s in samples if "gNEG4"   in s]
print(f"\n  Group-level mean CLR:")
print(f"    HnrnpM (n={len(group_HM)}): {clr_df[group_HM].values.mean():.4f}")
print(f"    NEG4   (n={len(group_NEG)}): {clr_df[group_NEG].values.mean():.4f}")

# Per-circRNA: how many have mean HM > mean NEG ?
mean_HM  = clr_df[group_HM ].mean(axis=1)
mean_NEG = clr_df[group_NEG].mean(axis=1)
above = (mean_HM > mean_NEG).sum()
total = len(clr_df)
print(f"\n  Per-circRNA: mean_HM > mean_NEG in {above:,} / {total:,} ({100*above/total:.1f}%)")
print(f"\nWrote: {OUT}")
