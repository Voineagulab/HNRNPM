"""
FinalAnalysis item 08 — cross-reference circRNA BSJ endpoints with rMATS AS-event
splice-site coordinates, using the new CIRIquant universe (5,820) and the
ADDITIVE BSJ-DE / FSJ-DE sets as the publication primary axes.

Run twice against TWO rMATS scopes:
  - PSD5+8 full     (publication splicing primary; ~2,617 sig events)
  - PSD5-only       (early-KD sensitivity;          ~2,393 sig events)

Match definition: same chr (rMATS 'chr' prefix stripped) + same strand + same
gene (rMATS GeneID present in circRNA's comma-separated gene_id list) + at least
one circRNA endpoint exactly matches an exon boundary in the rMATS event.

Strata tested:
  - universe         (5,820 circRNAs)
  - additive BSJ-DE  (56 circRNAs, publication primary for the circRNA axis)
  - additive FSJ-DE  (98 circRNAs, negative-control biology — should NOT enrich
                      for HnrnpM-differential AS, mirroring the CLIP test)
  - primary BSJ-DE   (58 circRNAs, sensitivity comparator)

Fisher's exact (one-sided "set enriched for matches_any") tests:
  - additive BSJ-DE vs universe
  - additive FSJ-DE vs universe   (negative-control expectation)
  - primary BSJ-DE  vs universe   (sensitivity)
For each rMATS scope (all / sig-full / sig-PSD5-only).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS")
OUT_DIR  = ROOT / "08_circRNA_AS_overlap"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BSJ_MAT  = ROOT / "01_circRNA_DE" / "bsj_matrix_CIRIquant_HnrnpM_PSD58.tsv"
BSJ_ADD  = ROOT / "01_circRNA_DE" / "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"
BSJ_PRI  = ROOT / "01_circRNA_DE" / "limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv"
FSJ_ADD  = ROOT / "01_circRNA_DE" / "FSJ" / "limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"
RMATS_FULL = ROOT / "07_Splicing_rMATS" / "rmats_raw"
RMATS_PSD5 = ROOT / "07_Splicing_rMATS" / "rmats_PSD5_only_HnrnpM"

EVENT_TYPES = ["SE", "RI", "A3SS", "A5SS", "MXE"]
START_COLS = {
    "SE":   ["exonStart_0base", "upstreamES", "downstreamES"],
    "RI":   ["riExonStart_0base", "upstreamES", "downstreamES"],
    "A3SS": ["longExonStart_0base", "shortES", "flankingES"],
    "A5SS": ["longExonStart_0base", "shortES", "flankingES"],
    "MXE":  ["1stExonStart_0base", "2ndExonStart_0base", "upstreamES", "downstreamES"],
}
END_COLS = {
    "SE":   ["exonEnd", "upstreamEE", "downstreamEE"],
    "RI":   ["riExonEnd", "upstreamEE", "downstreamEE"],
    "A3SS": ["longExonEnd", "shortEE", "flankingEE"],
    "A5SS": ["longExonEnd", "shortEE", "flankingEE"],
    "MXE":  ["1stExonEnd", "2ndExonEnd", "upstreamEE", "downstreamEE"],
}
FDR_THR  = 0.05
DPSI_THR = 0.10

# ---- Load circRNAs ----
circ_full = pd.read_csv(BSJ_MAT, sep="\t")
circ_full["chr"]    = circ_full["chr"].astype(str)
circ_full["start"]  = circ_full["start"].astype(int)
circ_full["end"]    = circ_full["end"].astype(int)
circ_full["strand"] = circ_full["strand"].astype(str)
circ_full["host_gene_list"] = circ_full["gene_id"].apply(
    lambda s: [] if pd.isna(s) else [x.strip() for x in str(s).split(",") if x.strip()])

bsj_add = pd.read_csv(BSJ_ADD, sep="\t")
bsj_pri = pd.read_csv(BSJ_PRI, sep="\t")
fsj_add = pd.read_csv(FSJ_ADD, sep="\t")

set_BSJ_add = set(bsj_add.loc[bsj_add["adj.P.Val"]<0.05, "circRNA"])
set_BSJ_pri = set(bsj_pri.loc[bsj_pri["adj.P.Val"]<0.05, "circRNA"])
set_FSJ_add = set(fsj_add.loc[fsj_add["adj.P.Val"]<0.05, "circRNA"])
print(f"Universe: {len(circ_full)}")
print(f"  BSJ_additive: {len(set_BSJ_add)} ; BSJ_primary: {len(set_BSJ_pri)} ; FSJ_additive: {len(set_FSJ_add)}")

# ---- Load rMATS events ----
def load_rmats(et, base_dir):
    df = pd.read_csv(base_dir / f"{et}.MATS.JCEC.txt", sep="\t",
                     dtype=str, low_memory=False)
    df["GeneID"]   = df["GeneID"].astype(str).str.strip('"')
    df["chr_norm"] = df["chr"].astype(str).str.replace("^chr","", regex=True)
    df["strand_norm"] = df["strand"]
    coord_cols = list(set(START_COLS[et]) | set(END_COLS[et]))
    for c in coord_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=coord_cols)
    for c in coord_cols:
        df[c] = df[c].astype(int)
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df["IncLevelDifference"] = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
    df["event_type"] = et
    return df

def filter_sig(df):
    return df[(df["FDR"].notna()) & (df["FDR"] < FDR_THR) &
              (df["IncLevelDifference"].notna()) &
              (df["IncLevelDifference"].abs() >= DPSI_THR)].copy()

def build_index(events_dict):
    idx = {et: {} for et in EVENT_TYPES}
    for et, df in events_dict.items():
        for _, r in df.iterrows():
            key = (r["chr_norm"], r["strand_norm"], r["GeneID"])
            if key not in idx[et]:
                idx[et][key] = (set(), set())
            for c in START_COLS[et]:
                idx[et][key][0].add(int(r[c]))
            for c in END_COLS[et]:
                idx[et][key][1].add(int(r[c]))
    return idx

print("\nLoading rMATS PSD5+8 (full) ...")
rmats_full = {et: load_rmats(et, RMATS_FULL) for et in EVENT_TYPES}
print("Loading rMATS PSD5-only ...")
rmats_psd5 = {et: load_rmats(et, RMATS_PSD5) for et in EVENT_TYPES}

idx_all       = build_index(rmats_full)
idx_sig_full  = build_index({et: filter_sig(df) for et, df in rmats_full.items()})
idx_sig_psd5  = build_index({et: filter_sig(df) for et, df in rmats_psd5.items()})
for label, ix in [("all", idx_all), ("sig_full", idx_sig_full), ("sig_PSD5_only", idx_sig_psd5)]:
    cnt = sum(len(d) for d in ix.values())
    print(f"  {label}: {cnt} (chr,strand,gene) entries indexed across event types")

def match_circ(row, idx):
    out = {et: False for et in EVENT_TYPES}
    for et in EVENT_TYPES:
        for g in row["host_gene_list"]:
            key = (row["chr"], row["strand"], g)
            if key in idx[et]:
                S, E = idx[et][key]
                if row["start"] in S or row["end"] in E:
                    out[et] = True; break
    return out

print("\nMatching circRNAs to all-rMATS / sig-PSD5+8 / sig-PSD5-only ...")
m_all       = circ_full.apply(lambda r: pd.Series(match_circ(r, idx_all)),      axis=1)
m_sig_full  = circ_full.apply(lambda r: pd.Series(match_circ(r, idx_sig_full)), axis=1)
m_sig_psd5  = circ_full.apply(lambda r: pd.Series(match_circ(r, idx_sig_psd5)), axis=1)

m_all.columns      = [f"matches_{et}_all"          for et in EVENT_TYPES]
m_sig_full.columns = [f"matches_{et}_sig_full"     for et in EVENT_TYPES]
m_sig_psd5.columns = [f"matches_{et}_sig_PSD5"     for et in EVENT_TYPES]
m_all["matches_any_all"]              = m_all.any(axis=1)
m_sig_full["matches_any_sig_full"]    = m_sig_full.any(axis=1)
m_sig_psd5["matches_any_sig_PSD5"]    = m_sig_psd5.any(axis=1)

result = pd.concat([circ_full, m_all, m_sig_full, m_sig_psd5], axis=1)
result["BSJ_additive_DE"] = result["circRNA"].isin(set_BSJ_add)
result["BSJ_primary_DE"]  = result["circRNA"].isin(set_BSJ_pri)
result["FSJ_additive_DE"] = result["circRNA"].isin(set_FSJ_add)

# Front column order
front = ["circRNA","gene_symbol","chr","start","end","strand","gene_id",
         "BSJ_additive_DE","BSJ_primary_DE","FSJ_additive_DE"]
match_cols = (["matches_any_all"] + [f"matches_{et}_all" for et in EVENT_TYPES]
              + ["matches_any_sig_full"] + [f"matches_{et}_sig_full" for et in EVENT_TYPES]
              + ["matches_any_sig_PSD5"] + [f"matches_{et}_sig_PSD5" for et in EVENT_TYPES])
remaining = [c for c in result.columns if c not in front + match_cols and c != "host_gene_list"]
result = result[front + match_cols + remaining]
result.to_csv(OUT_DIR / "circRNA_AS_overlap_PSD58.tsv", sep="\t", index=False)
print(f"\nWrote per-circRNA table: {OUT_DIR/'circRNA_AS_overlap_PSD58.tsv'}")

# ---- Summary per stratum × scope ----
def summary_for(df, stratum_name):
    rows = []
    for scope in ["all", "sig_full", "sig_PSD5"]:
        for et_label, col in [("any", f"matches_any_{scope}")] + \
                              [(et, f"matches_{et}_{scope}") for et in EVENT_TYPES]:
            n_total = len(df); n_match = int(df[col].sum())
            rows.append({
                "stratum": stratum_name, "rmats_scope": scope,
                "event_type": et_label,
                "n_total": n_total, "n_match": n_match,
                "pct_match": round(100*n_match/max(n_total,1), 2)
            })
    return rows

rows = []
rows += summary_for(result,                              "universe")
rows += summary_for(result[result["BSJ_additive_DE"]],   "BSJ_additive_DE")
rows += summary_for(result[result["BSJ_primary_DE"]],    "BSJ_primary_DE")
rows += summary_for(result[result["FSJ_additive_DE"]],   "FSJ_additive_DE")
summary = pd.DataFrame(rows)
summary.to_csv(OUT_DIR / "circRNA_AS_overlap_summary_PSD58.tsv", sep="\t", index=False)

# Pretty print: focus on the sig scopes
print("\n=== % matching at least one rMATS event of given event type ===")
for scope in ["all","sig_full","sig_PSD5"]:
    print(f"\n--- scope: {scope} ---")
    sub = summary[summary["rmats_scope"] == scope].pivot(
        index="event_type", columns="stratum", values="pct_match"
    ).reindex(["any","SE","RI","A3SS","A5SS","MXE"])[
        ["universe","BSJ_additive_DE","BSJ_primary_DE","FSJ_additive_DE"]]
    print(sub.to_string())

# ---- Fisher's exact ----
def fisher_for(scope, sub_df, bg_df, sub_name, bg_name):
    col = f"matches_any_{scope}"
    sub_match = sub_df[col].sum();  sub_n = len(sub_df)
    bg_match  = bg_df[col].sum();   bg_n  = len(bg_df)
    table = np.array([
        [sub_match, sub_n - sub_match],
        [bg_match - sub_match, (bg_n - sub_n) - (bg_match - sub_match)]
    ])
    or_, p_one = fisher_exact(table, alternative="greater")
    p_two = fisher_exact(table, alternative="two-sided")[1]
    return {"scope": scope, "comparison": f"{sub_name} vs {bg_name}",
            "sub_match": int(sub_match), "sub_n": sub_n,
            "bg_match": int(bg_match) - int(sub_match),
            "bg_n": bg_n - sub_n,
            "OR": round(or_, 3), "p_one_sided": p_one, "p_two_sided": p_two}

fisher_rows = []
for scope in ["all","sig_full","sig_PSD5"]:
    fisher_rows.append(fisher_for(scope, result[result["BSJ_additive_DE"]], result, "BSJ_additive_DE", "universe"))
    fisher_rows.append(fisher_for(scope, result[result["FSJ_additive_DE"]], result, "FSJ_additive_DE", "universe"))
    fisher_rows.append(fisher_for(scope, result[result["BSJ_primary_DE"]],  result, "BSJ_primary_DE",  "universe"))
fisher = pd.DataFrame(fisher_rows)
fisher.to_csv(OUT_DIR / "circRNA_AS_overlap_fisher_PSD58.tsv", sep="\t", index=False)
print("\n=== Fisher's exact (one-sided: enriched for matches_any vs universe) ===")
print(fisher.to_string(index=False))
print("\nDone.")
