"""
FinalAnalysis item 09 step 3 — Wilcoxon + Fisher tests for HnrnpM CLIP
flanking-intron enrichment, anchored on the ADDITIVE BSJ-DE set as the
publication primary. Also reports sensitivity numbers against the primary
~ group BSJ-DE set and against the additive FSJ-DE set.

Tests (per DE-defining axis):
  1. DE flank vs non-DE flank          (Wilcoxon rank-sum, one-sided greater)
  2. DE flank vs DE within-gene bg     (paired Wilcoxon, one-sided greater)
  3. DE flank vs non-host-gene bg      (Wilcoxon rank-sum, one-sided greater)
  4. Fisher's exact on bound_any_pos / bound_any_2x in DE-set vs rest of universe

Also writes the per-circRNA flank-score table with bound flags.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu, wilcoxon, fisher_exact

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV")
OUT_DIR = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "09_HnrnpM_CLIP_enrichment"
SIG = OUT_DIR / "clip_signal_per_region.tsv"

BSJ_ADD = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "01_circRNA_DE" / "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"
BSJ_PRI = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "01_circRNA_DE" / "limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv"
FSJ_ADD = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "01_circRNA_DE" / "FSJ" / "limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"

s = pd.read_csv(SIG, sep="\t")
bsj_add = pd.read_csv(BSJ_ADD, sep="\t")
bsj_pri = pd.read_csv(BSJ_PRI, sep="\t")
fsj_add = pd.read_csv(FSJ_ADD, sep="\t")

de_sets = {
    "BSJ_additive":   set(bsj_add.loc[bsj_add["adj.P.Val"] < 0.05, "circRNA"]),
    "BSJ_primary":    set(bsj_pri.loc[bsj_pri["adj.P.Val"] < 0.05, "circRNA"]),
    "FSJ_additive":   set(fsj_add.loc[fsj_add["adj.P.Val"] < 0.05, "circRNA"]),
}
for k, v in de_sets.items():
    print(f"DE set [{k}]: {len(v)} circRNAs")

flanks  = s[s["region_type"].isin(["flank5","flank3"])].copy()
within  = s[s["region_type"] == "within_gene_bg"].copy()
nonhost = s[s["region_type"] == "non_host_bg"].copy()

# Per-circRNA mean of flank5/flank3 + per-circRNA mean of within-gene bg
per_circ_flank  = flanks.groupby("circRNA")["log2_IP_over_Input"].mean().reset_index().rename(
    columns={"log2_IP_over_Input":"flank_score"})
per_circ_within = within.groupby("circRNA")["log2_IP_over_Input"].mean().reset_index().rename(
    columns={"log2_IP_over_Input":"within_gene_bg_score"})
m = per_circ_flank.merge(per_circ_within, on="circRNA", how="left")
print(f"\nPer-circRNA scores: {len(m)} with flank data; "
      f"{m['within_gene_bg_score'].notna().sum()} also with within-gene bg")
nonhost_score = nonhost["log2_IP_over_Input"].dropna().values

# ---- Wilcoxon tests for each DE-defining axis ----
results = []
for name, de_set in de_sets.items():
    m_local = m.copy()
    m_local["is_DE"] = m_local["circRNA"].isin(de_set)
    de_flank   = m_local.loc[ m_local["is_DE"], "flank_score"].dropna().values
    nde_flank  = m_local.loc[~m_local["is_DE"], "flank_score"].dropna().values
    if len(de_flank) < 3: continue
    # 1) DE flank vs non-DE flank
    u, p1 = mannwhitneyu(de_flank, nde_flank, alternative="greater")
    p1_two = mannwhitneyu(de_flank, nde_flank, alternative="two-sided")[1]
    results.append({
        "DE_axis": name, "test": "1. DE flank vs non-DE flank",
        "n1": len(de_flank), "n2": len(nde_flank),
        "median1": float(np.median(de_flank)), "median2": float(np.median(nde_flank)),
        "diff_medians": float(np.median(de_flank) - np.median(nde_flank)),
        "stat": float(u), "p_one_sided": p1, "p_two_sided": p1_two
    })
    # 2) DE flank vs DE within-gene bg (paired)
    de_pair = m_local[m_local["is_DE"] & m_local["within_gene_bg_score"].notna()].dropna(subset=["flank_score"])
    if len(de_pair) >= 5:
        w, p2 = wilcoxon(de_pair["flank_score"], de_pair["within_gene_bg_score"], alternative="greater")
        p2_two = wilcoxon(de_pair["flank_score"], de_pair["within_gene_bg_score"], alternative="two-sided")[1]
        results.append({
            "DE_axis": name, "test": "2. DE flank vs DE within-gene bg (paired)",
            "n1": len(de_pair), "n2": len(de_pair),
            "median1": float(np.median(de_pair["flank_score"])),
            "median2": float(np.median(de_pair["within_gene_bg_score"])),
            "diff_medians": float(np.median(de_pair["flank_score"]) - np.median(de_pair["within_gene_bg_score"])),
            "stat": float(w), "p_one_sided": p2, "p_two_sided": p2_two
        })
    # 3) DE flank vs non-host bg
    u3, p3 = mannwhitneyu(de_flank, nonhost_score, alternative="greater")
    p3_two = mannwhitneyu(de_flank, nonhost_score, alternative="two-sided")[1]
    results.append({
        "DE_axis": name, "test": "3. DE flank vs non-host-gene bg",
        "n1": len(de_flank), "n2": len(nonhost_score),
        "median1": float(np.median(de_flank)),
        "median2": float(np.median(nonhost_score)),
        "diff_medians": float(np.median(de_flank) - np.median(nonhost_score)),
        "stat": float(u3), "p_one_sided": p3, "p_two_sided": p3_two
    })

res = pd.DataFrame(results)
res.to_csv(OUT_DIR / "enrichment_test_summary.tsv", sep="\t", index=False, float_format="%.4g")
print("\n=== Wilcoxon tests ===")
for _, r in res.iterrows():
    print(f"  [{r['DE_axis']}] {r['test']}")
    print(f"    n1={r['n1']}, n2={r['n2']}, median {r['median1']:+.3f} vs {r['median2']:+.3f} (diff {r['diff_medians']:+.3f})")
    print(f"    P_one-sided = {r['p_one_sided']:.3g}; P_two-sided = {r['p_two_sided']:.3g}")

# ---- Per-circRNA bound flags + Fisher's exact ----
flank_wide = (flanks.pivot_table(index="circRNA", columns="region_type",
                                  values="log2_IP_over_Input", aggfunc="first")
                    .reset_index()
                    .rename(columns={"flank5":"flank5_log2","flank3":"flank3_log2"}))
flank_wide["max_flank_log2"]  = flank_wide[["flank5_log2","flank3_log2"]].max(axis=1)
flank_wide["mean_flank_log2"] = flank_wide[["flank5_log2","flank3_log2"]].mean(axis=1)
flank_wide["bound_any_pos"]  = flank_wide["max_flank_log2"]  > 0
flank_wide["bound_any_2x"]   = flank_wide["max_flank_log2"]  > 1
flank_wide["bound_mean_pos"] = flank_wide["mean_flank_log2"] > 0
flank_wide["bound_mean_2x"]  = flank_wide["mean_flank_log2"] > 1

# Annotate with additive BSJ-DE status (publication primary)
flank_wide["BSJ_add_logFC"]   = flank_wide["circRNA"].map(bsj_add.set_index("circRNA")["logFC"])
flank_wide["BSJ_add_adj_P"]   = flank_wide["circRNA"].map(bsj_add.set_index("circRNA")["adj.P.Val"])
flank_wide["BSJ_add_DE"]      = flank_wide["BSJ_add_adj_P"] < 0.05
flank_wide.to_csv(OUT_DIR / "circRNA_flank_score_PSD58.tsv", sep="\t", index=False, float_format="%.4g")

print("\n=== Fisher's exact for bound-flag enrichment ===")
fisher_rows = []
for name, de_set in de_sets.items():
    in_universe = flank_wide[flank_wide["bound_any_pos"].notna()]
    in_set      = in_universe[in_universe["circRNA"].isin(de_set)]
    rest        = in_universe[~in_universe["circRNA"].isin(de_set)]
    for flag in ["bound_any_pos","bound_any_2x"]:
        k = int(in_set[flag].sum()); n_set = len(in_set)
        K = int(rest[flag].sum());   n_rest = len(rest)
        if n_set == 0: continue
        or_val, p = fisher_exact([[k, n_set-k], [K, n_rest-K]], alternative="greater")
        fisher_rows.append({"DE_axis": name, "flag": flag,
                            "n_in_set": n_set, "n_universe_rest": n_rest,
                            "k_bound_in_set": k, "K_bound_rest": K,
                            "pct_bound_set":  round(100*k/n_set, 1),
                            "pct_bound_rest": round(100*K/n_rest, 1),
                            "OR": round(or_val, 3), "P_one_sided": p})
        print(f"  [{name}] {flag}: {k}/{n_set} = {100*k/n_set:.1f}% (set) vs {K}/{n_rest} = {100*K/n_rest:.1f}% (rest); OR={or_val:.2f}, P={p:.3g}")
fr = pd.DataFrame(fisher_rows)
fr.to_csv(OUT_DIR / "CLIP_bound_flag_fisher.tsv", sep="\t", index=False, float_format="%.4g")

print(f"\nWrote outputs to: {OUT_DIR}")
