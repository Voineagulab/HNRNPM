"""
FinalAnalysis item 06: compare ADDITIVE (~ group + PSD) gene-level DE to
Ho et al. 2021 (eLife 59654, fig4-data2). Publication primary is the additive
model — anchor the comparison there. Also runs the primary `~ group` set as
a sensitivity check.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr, chi2

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV")
HO_XLSX = ROOT / "PublicData" / "Ho_Elife_2021" / "elife-59654-fig4-data2-v2.xlsx"
LIMMA_ADD = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "02_GeneTx_DE" / "limma_gene_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"
LIMMA_PRI = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "02_GeneTx_DE" / "limma_gene_HnrnpM_vs_NEG_PSD58.tsv"
OUT_DIR   = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "06_Ho2021_Comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)
ALPHA = 0.05

# Ho fig4-data2
ho = pd.read_excel(HO_XLSX, sheet_name="Source Data 4", header=2).rename(columns={
    "Ensembl_gene_id":          "gene_id",
    "B7/SCR Log2(fold_change)": "Ho_B7_log2FC",
    "B7/SCR p_value":           "Ho_B7_p",
    "B7/SCR q_value":           "Ho_B7_q",
    "B7/SCR significant":       "Ho_B7_sig_str",
    "B9/SCR Log2(fold_change)": "Ho_B9_log2FC",
    "B9/SCR p_value":           "Ho_B9_p",
    "B9/SCR q_value":           "Ho_B9_q",
    "B9/SCR significant":       "Ho_B9_sig_str",
    "Bound":                    "Ho_Bound",
}).dropna(subset=["gene_id"])
ho["Ho_B7_sig"] = ho["Ho_B7_sig_str"].astype(str).str.lower() == "yes"
ho["Ho_B9_sig"] = ho["Ho_B9_sig_str"].astype(str).str.lower() == "yes"
ho["Ho_mean_log2FC"] = ho[["Ho_B7_log2FC","Ho_B9_log2FC"]].mean(axis=1)
ho["Ho_sig_strict"] = (ho["Ho_B7_sig"] & ho["Ho_B9_sig"] &
    (np.sign(ho["Ho_B7_log2FC"]) == np.sign(ho["Ho_B9_log2FC"])) &
    (ho["Ho_B7_log2FC"] != 0))
print(f"Ho rows: {len(ho)} | strict-sig: {ho['Ho_sig_strict'].sum()}")

def run(limma_path, label, out_per_gene, out_summary):
    lm = pd.read_csv(limma_path, sep="\t").rename(columns={
        "logFC": "limma_logFC", "adj.P.Val": "limma_adj_P"
    })[["gene_id","gene_name","limma_logFC","limma_adj_P"]]
    m = ho[["gene_id","Ho_B7_log2FC","Ho_B9_log2FC","Ho_mean_log2FC",
            "Ho_B7_q","Ho_B9_q","Ho_sig_strict","Ho_Bound"]].merge(lm, on="gene_id", how="outer")
    m["Ho_sig_strict"] = m["Ho_sig_strict"].fillna(False)
    m["limma_sig"]     = (m["limma_adj_P"] < ALPHA).fillna(False)

    def direction(ho_lfc, our_lfc, hs, os_):
        if not (hs and os_): return "NS"
        if pd.isna(ho_lfc) or pd.isna(our_lfc): return "NS"
        if ho_lfc > 0 and our_lfc > 0: return "concordant_UP"
        if ho_lfc < 0 and our_lfc < 0: return "concordant_DOWN"
        return "discordant"
    m["direction_limma_vs_Ho"] = [
        direction(h, o, hs, os_) for h, o, hs, os_ in
        zip(m["Ho_mean_log2FC"], m["limma_logFC"], m["Ho_sig_strict"], m["limma_sig"])]
    m = m.sort_values(["Ho_sig_strict","limma_adj_P"], ascending=[False, True])
    m.to_csv(out_per_gene, sep="\t", index=False)

    universe = m[m["limma_adj_P"].notna()]
    N = len(universe)
    ho_in_uni = universe[universe["Ho_sig_strict"]]
    K_up   = int((ho_in_uni["Ho_mean_log2FC"] > 0).sum())
    K_down = int((ho_in_uni["Ho_mean_log2FC"] < 0).sum())
    our_sig = universe[universe["limma_adj_P"] < ALPHA]
    n_up   = int((our_sig["limma_logFC"] > 0).sum())
    n_down = int((our_sig["limma_logFC"] < 0).sum())
    overlap = our_sig[our_sig["Ho_sig_strict"]]
    k_up   = int(((overlap["Ho_mean_log2FC"] > 0) & (overlap["limma_logFC"] > 0)).sum())
    k_down = int(((overlap["Ho_mean_log2FC"] < 0) & (overlap["limma_logFC"] < 0)).sum())
    k_disc = int(((np.sign(overlap["Ho_mean_log2FC"]) != np.sign(overlap["limma_logFC"])) &
                  (overlap["Ho_mean_log2FC"] != 0) & (overlap["limma_logFC"] != 0)).sum())

    P_up   = hypergeom.sf(k_up - 1,   N, K_up,   n_up)
    P_down = hypergeom.sf(k_down - 1, N, K_down, n_down)
    P_combined = chi2.sf(-2 * (np.log(P_up) + np.log(P_down)), df = 4)
    fold_up   = (k_up   / max(n_up, 1))   / (max(K_up, 1)   / max(N,1))
    fold_down = (k_down / max(n_down, 1)) / (max(K_down, 1) / max(N,1))
    rho_K, p_K = spearmanr(ho_in_uni["Ho_mean_log2FC"], ho_in_uni["limma_logFC"], nan_policy="omit")
    rho_k, p_k = spearmanr(overlap["Ho_mean_log2FC"], overlap["limma_logFC"], nan_policy="omit") \
                 if len(overlap) >= 3 else (np.nan, np.nan)

    stats = pd.DataFrame([{
        "model": label, "N": N,
        "K_up": K_up, "K_down": K_down,
        "n_up": n_up, "n_down": n_down,
        "k_up": k_up, "k_down": k_down, "k_disc": k_disc,
        "fold_enrichment_UP": round(fold_up, 2), "fold_enrichment_DOWN": round(fold_down, 2),
        "P_up": P_up, "P_down": P_down, "P_combined": P_combined,
        "Spearman_rho_K_genes": round(rho_K, 3), "Spearman_p_K_genes": p_K,
        "Spearman_rho_overlap": round(rho_k, 3) if not np.isnan(rho_k) else np.nan,
        "Spearman_p_overlap":   p_k,
    }])
    stats.to_csv(out_summary, sep="\t", index=False)
    print(f"\n=== [{label}] direction-stratified hypergeometric ===")
    print(f"  N = {N}; K (UP/DOWN) = {K_up}/{K_down}; n (UP/DOWN) = {n_up}/{n_down}")
    print(f"  k (concordant UP/DOWN) = {k_up}/{k_down}; k_disc = {k_disc}")
    print(f"  Fold enrichment: UP {fold_up:.2f}, DOWN {fold_down:.2f}")
    print(f"  P_combined (Fisher's method) = {P_combined:.2e}")
    print(f"  Spearman rho on overlap = {rho_k:.3f}")
    return stats

# Publication primary: additive
add_stats = run(LIMMA_ADD, "additive ~group+PSD (publication primary)",
                OUT_DIR / "Ho2021_vs_ours_HnrnpM_DE_PSD58.tsv",
                OUT_DIR / "Ho2021_vs_ours_HnrnpM_summary_PSD58.tsv")
# Sensitivity: primary ~group
pri_stats = run(LIMMA_PRI, "primary ~group (sensitivity)",
                OUT_DIR / "Ho2021_vs_ours_HnrnpM_DE_PSD58_primaryGroup.tsv",
                OUT_DIR / "Ho2021_vs_ours_HnrnpM_summary_PSD58_primaryGroup.tsv")

combined = pd.concat([add_stats, pri_stats], ignore_index=True)
combined.to_csv(OUT_DIR / "Ho2021_vs_ours_HnrnpM_summary_both_models.tsv",
                sep="\t", index=False, float_format="%.4g")
print("\nWrote outputs to:", OUT_DIR)
