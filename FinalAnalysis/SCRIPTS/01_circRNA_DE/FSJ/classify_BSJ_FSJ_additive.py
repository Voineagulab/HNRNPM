"""
Cross-classify BSJ-DE circRNAs against their FSJ status — sensitivity version
using the additive (`~ group + PSD`) model on both sides.

Compares to the primary (`~ group`) classification to assess whether the
86% "independent back-splicing" headline holds up when PSD is included as
covariate in both DE tests.
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE")
OUT  = ROOT / "FSJ"

# Inputs
bsj_pri  = pd.read_csv(ROOT / "limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv",         sep="\t")
bsj_add  = pd.read_csv(ROOT / "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv",  sep="\t")
fsj_pri  = pd.read_csv(OUT  / "limma_FSJ_HnrnpM_vs_NEG_PSD58.tsv",         sep="\t")
fsj_add  = pd.read_csv(OUT  / "limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv",  sep="\t")
fsj_fail = pd.read_csv(OUT  / "fsj_filter_fail_circRNAs.tsv",              sep="\t")
fsj_fail_set = set(fsj_fail["circRNA"])

LFC_T_SOFT = np.log2(1.5)   # |logFC| threshold for the soft FSJ-sig rule

def classify(bsj_df, fsj_df, label):
    # BSJ-sig: adj.P < 0.05 (strict, as in main DE)
    # FSJ-sig: SOFT rule — raw P < 0.05 AND |logFC_FSJ| >= log2(1.5).
    # Rationale: the per-axis FSJ FDR test is conservative for low-coverage
    # events (FSJ AveExpr ~ -5 for many circRNAs flanking BSJ-DE loci), so
    # genuine large-FC FSJ changes don't survive BH at ~5,550 tests.
    # We use the raw-P + magnitude criterion to capture biologically meaningful
    # FSJ shifts while still requiring statistical evidence.
    bsj_keep = bsj_df[["circRNA","gene_symbol","logFC","adj.P.Val"]].rename(
        columns={"logFC":"logFC_BSJ","adj.P.Val":"adjP_BSJ"})
    fsj_keep = fsj_df[["circRNA","logFC","P.Value","adj.P.Val"]].rename(
        columns={"logFC":"logFC_FSJ","P.Value":"rawP_FSJ","adj.P.Val":"adjP_FSJ"})
    m = bsj_keep.merge(fsj_keep, on="circRNA", how="left")
    m["sig_BSJ"] = m["adjP_BSJ"] < 0.05
    m["sig_FSJ"] = (m["rawP_FSJ"] < 0.05) & (m["logFC_FSJ"].abs() >= LFC_T_SOFT)
    m["dir_BSJ"] = np.sign(m["logFC_BSJ"])
    m["dir_FSJ"] = np.sign(m["logFC_FSJ"])
    m["in_fsj_test"] = m["adjP_FSJ"].notna()
    m["fsj_filter_fail"] = m["circRNA"].isin(fsj_fail_set)

    def cls(row):
        if row["sig_BSJ"] and row["fsj_filter_fail"]:
            return "extreme_linear_collapse"
        if row["sig_BSJ"] and not row["in_fsj_test"]:
            return "extreme_linear_collapse"
        if row["sig_BSJ"] and row["sig_FSJ"]:
            if row["dir_BSJ"] == row["dir_FSJ"]:
                return "co_regulated"
            else:
                return "opposite_regulation"
        if row["sig_BSJ"] and not row["sig_FSJ"]:
            return "independent_backsplicing"
        if not row["sig_BSJ"] and row["sig_FSJ"]:
            return "FSJ_only_NS_BSJ"
        return "neither"
    m["class"] = m.apply(cls, axis=1)

    bsj_sig = m[m["sig_BSJ"]]
    counts = bsj_sig["class"].value_counts().reindex(
        ["independent_backsplicing","co_regulated","opposite_regulation","extreme_linear_collapse"],
        fill_value=0)
    print(f"\n=== {label} model — BSJ-DE circRNA classification ===")
    print(f"  Total BSJ-DE: {len(bsj_sig)}  ({m['sig_BSJ'].sum()} primary-or-additive)")
    for k, v in counts.items():
        pct = 100*v/len(bsj_sig) if len(bsj_sig) > 0 else 0
        print(f"  {k:<28} {v:>3}  ({pct:.1f}%)")
    return m, counts

m_pri, counts_pri = classify(bsj_pri, fsj_pri, "PRIMARY ~group")
m_add, counts_add = classify(bsj_add, fsj_add, "ADDITIVE ~group + PSD")

# Side-by-side summary
summary = pd.DataFrame({
    "primary": counts_pri,
    "additive": counts_add,
}).reset_index().rename(columns={"index":"class"})
summary.to_csv(OUT / "BSJ_FSJ_classification_primary_vs_additive_summary.tsv", sep="\t", index=False)

# Save the additive-model classification per circRNA
m_add_out = m_add[["circRNA","gene_symbol","logFC_BSJ","adjP_BSJ","logFC_FSJ","rawP_FSJ","adjP_FSJ",
                    "sig_BSJ","sig_FSJ","in_fsj_test","fsj_filter_fail","class"]].copy()
m_add_out.to_csv(OUT / "BSJ_FSJ_classification_additive.tsv", sep="\t", index=False)

# How many circRNAs change class between primary and additive?
m_compare = m_pri[["circRNA","class"]].rename(columns={"class":"class_primary"}).merge(
            m_add[["circRNA","class"]].rename(columns={"class":"class_additive"}),
            on="circRNA")
m_compare_bsj_sig = m_compare[m_compare["circRNA"].isin(set(bsj_pri.loc[bsj_pri["adj.P.Val"]<0.05,"circRNA"]) |
                                                          set(bsj_add.loc[bsj_add["adj.P.Val"]<0.05,"circRNA"]))]
class_changes = m_compare_bsj_sig[m_compare_bsj_sig["class_primary"] != m_compare_bsj_sig["class_additive"]]
print(f"\n--- Per-circRNA class changes between primary and additive (BSJ-sig in either) ---")
print(f"Total compared: {len(m_compare_bsj_sig)}; class changes: {len(class_changes)}")
if len(class_changes) > 0:
    print(class_changes.to_string(index=False))
    class_changes.to_csv(OUT / "BSJ_FSJ_classification_primary_vs_additive_diffs.tsv",
                          sep="\t", index=False)

# Print final summary table
print(f"\n=== Final summary table ===")
print(summary.to_string(index=False))
