"""
FinalAnalysis item 06b — bound vs unbound replication test using the ADDITIVE
gene-level DE (publication primary).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/06_Ho2021_Comparison")
PER_GENE = ROOT / "Ho2021_vs_ours_HnrnpM_DE_PSD58.tsv"
OUT_TSV  = ROOT / "bound_vs_unbound_replication_PSD58.tsv"
ALPHA = 0.05

m = pd.read_csv(PER_GENE, sep="\t", low_memory=False)
def to_bool(s):
    return s.map({True: True, False: False, "True": True, "False": False})
m["Ho_sig_strict"] = to_bool(m["Ho_sig_strict"])
m["Ho_Bound"]      = to_bool(m["Ho_Bound"])
m["limma_sig"]     = to_bool(m["limma_sig"])

ho = m[m["Ho_sig_strict"].fillna(False) & m["Ho_Bound"].notna()].copy()
print(f"Ho strict-sig with non-NA Bound: {len(ho)}")
print(ho["Ho_Bound"].value_counts())

sub = ho[ho["limma_adj_P"].notna()].copy()
sub["our_sig"] = sub["limma_sig"].fillna(False)
sub["our_concordant"] = (
    sub["our_sig"] &
    (np.sign(sub["limma_logFC"]) == np.sign(sub["Ho_mean_log2FC"])) &
    (sub["limma_logFC"] != 0) & (sub["Ho_mean_log2FC"] != 0)
)
bound   = sub[sub["Ho_Bound"] == True]
unbound = sub[sub["Ho_Bound"] == False]

table = np.array([
    [bound["our_concordant"].sum(),   (~bound["our_concordant"]).sum()],
    [unbound["our_concordant"].sum(), (~unbound["our_concordant"]).sum()],
])
or_c, p_c_one = fisher_exact(table, alternative="greater")
p_c_two = fisher_exact(table, alternative="two-sided")[1]

table_sig = np.array([
    [bound["our_sig"].sum(),   (~bound["our_sig"]).sum()],
    [unbound["our_sig"].sum(), (~unbound["our_sig"]).sum()],
])
or_s, p_s_one = fisher_exact(table_sig, alternative="greater")
p_s_two = fisher_exact(table_sig, alternative="two-sided")[1]

rho_b, p_rho_b = spearmanr(bound["Ho_mean_log2FC"],   bound["limma_logFC"],   nan_policy="omit") \
                 if len(bound)   >= 3 else (np.nan, np.nan)
rho_u, p_rho_u = spearmanr(unbound["Ho_mean_log2FC"], unbound["limma_logFC"], nan_policy="omit") \
                 if len(unbound) >= 3 else (np.nan, np.nan)

rows = []
for label, s in [("Bound", bound), ("Unbound", unbound)]:
    rows.append({
        "model":"additive (~group+PSD)", "stratum": label,
        "n_total": len(s), "n_our_sig": int(s["our_sig"].sum()),
        "pct_our_sig": round(100 * s["our_sig"].sum() / max(len(s), 1), 2),
        "n_concordant": int(s["our_concordant"].sum()),
        "pct_concordant": round(100 * s["our_concordant"].sum() / max(len(s), 1), 2),
        "n_discordant": int((s["our_sig"] & ~s["our_concordant"]).sum()),
    })
pd.DataFrame(rows).to_csv(OUT_TSV, sep="\t", index=False)

print("\n=== Per-stratum counts (additive) ===")
print(pd.DataFrame(rows).to_string(index=False))
print(f"\n=== Bound vs Unbound (Fisher one-sided) ===")
print(f"  Concordant: bound {bound['our_concordant'].sum()}/{len(bound)} ({100*bound['our_concordant'].sum()/len(bound):.1f}%) vs "
      f"unbound {unbound['our_concordant'].sum()}/{len(unbound)} ({100*unbound['our_concordant'].sum()/len(unbound):.1f}%); "
      f"OR={or_c:.3f}, P_one={p_c_one:.4f}, P_two={p_c_two:.4f}")
print(f"  Any-sig:    bound {bound['our_sig'].sum()}/{len(bound)} vs unbound {unbound['our_sig'].sum()}/{len(unbound)}; "
      f"OR={or_s:.3f}, P_one={p_s_one:.4f}")
print(f"  Spearman rho (bound)   = {rho_b:.3f} (p={p_rho_b:.2e})")
print(f"  Spearman rho (unbound) = {rho_u:.3f} (p={p_rho_u:.2e})")
print(f"\nWrote: {OUT_TSV}")
