"""
FinalAnalysis item 09 step 2 — extract per-region mean CLIP coverage from the
4 hg38 HnrnpM eCLIP BigWigs (2 IP + 2 Input) over flanking introns and
background regions. Library-size-normalise, then compute per-region
log2((IP+eps) / (Input+eps)) with eps = 5th percentile of non-zero signal.
"""
import numpy as np
import pandas as pd
from pathlib import Path
import pyBigWig

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV")
BW_DIR  = ROOT / "PublicData" / "Ho_Elife_2021" / "hg38"
OUT_DIR = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "09_HnrnpM_CLIP_enrichment"

BWS = {
    "IP1":    BW_DIR / "GSM3119363_Hnrnpm_rep1.hg38.bw",
    "IP2":    BW_DIR / "GSM3119364_Hnrnpm_rep2.hg38.bw",
    "Input1": BW_DIR / "GSM3119365_Input_rep1.hg38.bw",
    "Input2": BW_DIR / "GSM3119366_Input_rep2.hg38.bw",
}

bw_handles = {}; totals = {}
for k, p in BWS.items():
    bw = pyBigWig.open(str(p))
    h = bw.header()
    totals[k] = h.get("sumData", None)
    if totals[k] is None:
        totals[k] = sum(bw.stats(c, type="sum")[0] or 0 for c in bw.chroms())
    bw_handles[k] = bw
    print(f"{k}: total signal = {totals[k]:.3e}")

total_IP    = (totals["IP1"]    + totals["IP2"])    / 2
total_Input = (totals["Input1"] + totals["Input2"]) / 2

def norm_mean(bw, total, chrom, s, e):
    v = bw.stats(chrom, int(s), int(e), type="mean")[0]
    if v is None: return 0.0
    return v / total * 1e9

def query_region(chrom, s, e):
    cc = "chr" + str(chrom) if not str(chrom).startswith("chr") else str(chrom)
    if cc not in bw_handles["IP1"].chroms():
        return (np.nan, np.nan)
    ip1 = norm_mean(bw_handles["IP1"],    total_IP,    cc, s, e)
    ip2 = norm_mean(bw_handles["IP2"],    total_IP,    cc, s, e)
    in1 = norm_mean(bw_handles["Input1"], total_Input, cc, s, e)
    in2 = norm_mean(bw_handles["Input2"], total_Input, cc, s, e)
    return ((ip1 + ip2) / 2, (in1 + in2) / 2)

flank   = pd.read_csv(OUT_DIR / "flanking_introns_per_circRNA.tsv", sep="\t")
within  = pd.read_csv(OUT_DIR / "within_gene_background_regions.tsv", sep="\t")
nonhost = pd.read_csv(OUT_DIR / "non_host_gene_background_regions.tsv", sep="\t")
print(f"\nFlanks rows: {len(flank)}; within-gene bg: {len(within)}; non-host bg: {len(nonhost)}")

def query_per_circ_flanks(df):
    rows = []
    for _, r in df.iterrows():
        for which in ["flank5","flank3"]:
            s = r[f"{which}_start"]; e = r[f"{which}_end"]
            if pd.isna(s) or pd.isna(e): continue
            ip, inp = query_region(r["chr"], int(s), int(e))
            rows.append({"circRNA": r["circRNA"], "gene_symbol": r["gene_symbol"],
                         "region_type": which, "chr": r["chr"],
                         "start": int(s), "end": int(e), "length": int(e-s),
                         "IP_signal": ip, "Input_signal": inp})
    return pd.DataFrame(rows)

print("Querying flanking introns ...")
flank_sig = query_per_circ_flanks(flank); print(f"  flank rows: {len(flank_sig)}")

print("Querying within-gene background regions ...")
w_rows = []
for _, r in within.iterrows():
    ip, inp = query_region(r["chr"], int(r["start"]), int(r["end"]))
    w_rows.append({"circRNA": r["circRNA"], "region_type": "within_gene_bg",
                   "chr": r["chr"], "start": int(r["start"]), "end": int(r["end"]),
                   "length": int(r["end"]-r["start"]),
                   "IP_signal": ip, "Input_signal": inp})
within_sig = pd.DataFrame(w_rows)

print("Querying non-host-gene background regions ...")
n_rows = []
for _, r in nonhost.iterrows():
    ip, inp = query_region(r["chr"], int(r["start"]), int(r["end"]))
    n_rows.append({"circRNA": np.nan, "region_type": "non_host_bg",
                   "chr": r["chr"], "start": int(r["start"]), "end": int(r["end"]),
                   "length": int(r["end"]-r["start"]),
                   "IP_signal": ip, "Input_signal": inp})
nonhost_sig = pd.DataFrame(n_rows)

all_sig = pd.concat([flank_sig, within_sig, nonhost_sig], ignore_index=True)
all_sig = all_sig.dropna(subset=["IP_signal","Input_signal"])
nz = all_sig[(all_sig["IP_signal"]>0) & (all_sig["Input_signal"]>0)]
eps = float(np.percentile(nz[["IP_signal","Input_signal"]].values.flatten(), 5))
print(f"\nPseudocount eps: {eps:.4g}")

all_sig["log2_IP_over_Input"] = np.log2((all_sig["IP_signal"] + eps) / (all_sig["Input_signal"] + eps))

OUT_TSV = OUT_DIR / "clip_signal_per_region.tsv"
all_sig.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\nWrote: {OUT_TSV} ({len(all_sig)} rows)")

print("\nMedian log2(IP/Input) by region type:")
print(all_sig.groupby("region_type")["log2_IP_over_Input"].agg(["count","median","mean"]).to_string())

for k, bw in bw_handles.items(): bw.close()
