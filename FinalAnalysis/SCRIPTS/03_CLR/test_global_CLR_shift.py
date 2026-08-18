"""
Statistical tests for the global CLR shift between HnrnpM and NEG4.
Three complementary tests:
  1. Per-sample mean CLR: t-test + Wilcoxon on the across-circRNA average CLR
     per sample (4 NEG vs 3 HnrnpM).
  2. Per-circRNA paired test: Wilcoxon signed-rank + paired t-test on
     (mean_HM - mean_NEG) across the 5,820 circRNAs.
  3. Sign / binomial test on the fraction-above-diagonal (3,970 / 5,820).
  4. Chi-squared on the 4-bin CLR distribution (NEG vs HnrnpM).
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu, ttest_rel, wilcoxon, binomtest, chi2_contingency

CLR_F = "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/03_CLR/clr_matrix_HnrnpM_PSD58.tsv"

dat = pd.read_csv(CLR_F, sep="\t")
meta_cols = ["circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id"]
samples = [c for c in dat.columns if c not in meta_cols]
mat = dat[samples].values

HM  = ["5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22"]
NEG = ["9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5"]

# --- Test 1: per-sample mean CLR ---
sample_means = mat.mean(axis=0)
sample_df = pd.DataFrame({"sample": samples, "mean_CLR": sample_means,
                          "group": ["HnrnpM" if "gHnrnpM" in s else "NEG4" for s in samples]})
hm_means  = sample_df.loc[sample_df.group=="HnrnpM", "mean_CLR"].values
neg_means = sample_df.loc[sample_df.group=="NEG4",   "mean_CLR"].values

t_stat, t_p = ttest_ind(hm_means, neg_means, equal_var=False, alternative="greater")
u_stat, u_p = mannwhitneyu(hm_means, neg_means, alternative="greater")

print("=== Test 1: per-sample mean CLR (4 NEG vs 3 HnrnpM, across-circRNA average) ===")
print(sample_df.to_string(index=False))
print(f"  HnrnpM mean of means: {hm_means.mean():.4f}  (n=3)")
print(f"  NEG4   mean of means: {neg_means.mean():.4f}  (n=4)")
print(f"  Welch t-test (greater): t = {t_stat:.3f}, P = {t_p:.4g}")
print(f"  Mann-Whitney U (greater): U = {u_stat:.0f}, P = {u_p:.4g}")

# --- Test 2: per-circRNA paired test ---
m_HM  = mat[:, [samples.index(s) for s in HM ]].mean(axis=1)
m_NEG = mat[:, [samples.index(s) for s in NEG]].mean(axis=1)
diff = m_HM - m_NEG
# Drop ties for Wilcoxon (zero diffs)
nonzero_diff = diff[diff != 0]
w_stat, w_p = wilcoxon(nonzero_diff, alternative="greater")
tp_stat, tp_p = ttest_rel(m_HM, m_NEG, alternative="greater")

print("\n=== Test 2: per-circRNA paired test on (mean_HM - mean_NEG) ===")
print(f"  n circRNAs:                  {len(diff):,}")
print(f"  n with non-zero difference: {len(nonzero_diff):,}")
print(f"  Median difference:           {np.median(diff):+.5f}")
print(f"  Mean   difference:           {diff.mean():+.5f}")
print(f"  Wilcoxon signed-rank (greater): W = {w_stat:.3g}, P = {w_p:.4g}")
print(f"  Paired t-test (greater):       t = {tp_stat:.3f}, P = {tp_p:.4g}")

# --- Test 3: sign / binomial test on direction ---
n_above = int((m_HM > m_NEG).sum())
n_below = int((m_HM < m_NEG).sum())
n_tied  = int((m_HM == m_NEG).sum())
n_total = len(m_HM)
bt = binomtest(n_above, n_above + n_below, p=0.5, alternative="greater")
print("\n=== Test 3: sign / binomial test on direction (above y=x diagonal) ===")
print(f"  Above diagonal: {n_above:,} / {n_above+n_below:,} ({100*n_above/(n_above+n_below):.1f}%)")
print(f"  Below diagonal: {n_below:,} ; tied: {n_tied:,}")
print(f"  Binomial (one-sided, p>0.5): P = {bt.pvalue:.4g}")

# --- Test 4: Chi-squared on bin distribution ---
bins = [(0, 0.1), (0.1, 0.5), (0.5, 0.9), (0.9, 1.001)]
labels = ["<0.1","0.1-0.5","0.5-0.9",">0.9"]
def bin_counts(values):
    v = values[values > 0]  # detected only, same as the density plot
    out = []
    for lo, hi in bins:
        out.append(int(((v >= lo) & (v < hi)).sum()))
    return out
hm_vals  = mat[:, [samples.index(s) for s in HM]].flatten()
neg_vals = mat[:, [samples.index(s) for s in NEG]].flatten()
ct = np.array([bin_counts(neg_vals), bin_counts(hm_vals)])
chi2, p_chi, dof, exp = chi2_contingency(ct)
print("\n=== Test 4: Chi-squared on 4-bin distribution (detected entries) ===")
df_ct = pd.DataFrame(ct, columns=labels, index=["NEG4","HnrnpM"])
print(df_ct.to_string())
print(f"  Chi-squared = {chi2:.1f}, dof = {dof}, P = {p_chi:.4g}")

# --- Save summary ---
out = pd.DataFrame([
    {"test":"per_sample_mean_CLR  Welch_t (HM>NEG)",
     "stat":t_stat, "p":t_p, "df":"n1=3, n2=4"},
    {"test":"per_sample_mean_CLR  Mann-Whitney_U (HM>NEG)",
     "stat":u_stat, "p":u_p, "df":"n1=3, n2=4"},
    {"test":"per_circRNA_paired  Wilcoxon (HM>NEG)",
     "stat":float(w_stat), "p":float(w_p), "df":f"n={len(nonzero_diff)}"},
    {"test":"per_circRNA_paired  paired_t (HM>NEG)",
     "stat":tp_stat, "p":tp_p, "df":f"n={n_total}"},
    {"test":"sign_test  (above y=x)",
     "stat":n_above, "p":bt.pvalue,
     "df":f"trials={n_above+n_below}, k={n_above}"},
    {"test":"bin_distribution  Chi-squared",
     "stat":chi2, "p":p_chi, "df":f"dof={dof}"},
])
out.to_csv("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/03_CLR/global_CLR_shift_tests.tsv",
           sep="\t", index=False, float_format="%.4g")
print(f"\nWrote: global_CLR_shift_tests.tsv")
