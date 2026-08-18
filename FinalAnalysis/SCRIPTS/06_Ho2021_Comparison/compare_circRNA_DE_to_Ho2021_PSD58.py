"""
Compare our FinalAnalysis circRNA DE (additive ~ group + PSD primary) to Ho et al. 2021
circRNA-sig set (`elife-59654-fig5-data1-v2.xlsx`, Source Data 6).

Ho coordinates are hg19 — lift to hg38 with CrossMap before joining on coordinates.
Two comparison strategies, in order of stringency:

  (A) Gene-symbol overlap: which of our 56 additive BSJ-DE host genes appear in
      Ho's 264 unique host genes? Doesn't require coordinate matching, gives a
      "do the same loci come up?" answer.

  (B) Coordinate overlap on hg38: after lifting Ho's BSJ coords to hg38, match
      each Ho circRNA to our circRNA universe by exact (chr, start, end, strand).
      Reports overlaps in:
        - our full 5,820 universe
        - our additive BSJ-DE set (56)
        - our additive FSJ-DE set (98)

Hypergeometric test for the BSJ-DE × Ho-circRNA-DE intersection, using our
5,820 universe as N and Ho-circRNAs-mapped-into-it as K.
"""
import re
import subprocess
import tempfile
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, hypergeom

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV")
HO_XLSX = ROOT / "PublicData" / "Ho_Elife_2021" / "elife-59654-fig5-data1-v2.xlsx"
CHAIN   = "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/hg19ToHg38.over.chain"
BSJ_ADD = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "01_circRNA_DE" / "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"
FSJ_ADD = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "01_circRNA_DE" / "FSJ" / "limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"
BSJ_MAT = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "01_circRNA_DE" / "bsj_matrix_CIRIquant_HnrnpM_PSD58.tsv"
OUT_DIR = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "06_Ho2021_Comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- Load Ho circRNA list ----
ho = pd.read_excel(HO_XLSX, sheet_name="Source Data 6", header=2)
def parse(s):
    if not isinstance(s, str): return None
    m = re.match(r'^([^_]+)_(chr[^_]+)_(\d+)_(\d+)_([+-])$', s)
    return m.groups() if m else None
parsed = ho["circRNA"].apply(parse)
ho["gene_symbol"] = parsed.apply(lambda t: t[0] if t else None)
ho["hg19_chr"]    = parsed.apply(lambda t: t[1] if t else None)
ho["hg19_start"]  = parsed.apply(lambda t: int(t[2]) if t else None)
ho["hg19_end"]    = parsed.apply(lambda t: int(t[3]) if t else None)
ho["strand"]      = parsed.apply(lambda t: t[4] if t else None)
ho["Ho_mean_log2FC"] = ho[["log2FoldChange B7/SCR","log2FoldChange B9/SCR"]].mean(axis=1)
print(f"Ho circRNA rows: {len(ho)} | unique gene symbols: {ho['gene_symbol'].nunique()}")
print(f"  Bound (CLIP): {ho['bound'].sum()} ; unbound: {(~ho['bound']).sum()}")

# ---- Load our DE results + universe ----
bsj_add = pd.read_csv(BSJ_ADD, sep="\t")
fsj_add = pd.read_csv(FSJ_ADD, sep="\t")
universe = pd.read_csv(BSJ_MAT, sep="\t")
universe["coord_key"] = universe.apply(
    lambda r: f"{r['chr']}:{int(r['start'])}-{int(r['end'])}:{r['strand']}", axis=1)
print(f"\nOur universe: {len(universe)} circRNAs")
print(f"Our additive BSJ-DE: {(bsj_add['adj.P.Val']<0.05).sum()} circRNAs")
print(f"Our additive FSJ-DE: {(fsj_add['adj.P.Val']<0.05).sum()} circRNAs")

# ---- (A) Gene-symbol overlap ----
our_universe_genes = set(universe["gene_symbol"].dropna())
our_bsj_de_genes   = set(bsj_add.loc[bsj_add["adj.P.Val"]<0.05, "gene_symbol"].dropna())
our_fsj_de_genes   = set(fsj_add.loc[fsj_add["adj.P.Val"]<0.05, "gene_symbol"].dropna())
ho_genes           = set(ho["gene_symbol"].dropna())
sym_a = our_universe_genes & ho_genes
sym_b = our_bsj_de_genes   & ho_genes
sym_c = our_fsj_de_genes   & ho_genes
print(f"\n--- (A) Gene-symbol overlap ---")
print(f"Ho gene symbols in our universe:    {len(sym_a)} / {len(ho_genes)} = {100*len(sym_a)/len(ho_genes):.1f}%")
print(f"Ho symbols also in our BSJ-DE host genes: {len(sym_b)} (BSJ-DE host gene count: {len(our_bsj_de_genes)})")
print(f"Ho symbols also in our FSJ-DE host genes: {len(sym_c)} (FSJ-DE host gene count: {len(our_fsj_de_genes)})")
if sym_b: print(f"  Shared between BSJ-DE & Ho:\n    {sorted(sym_b)}")

# ---- (B) Coordinate lift hg19 -> hg38, then exact match ----
print("\n--- (B) Coordinate lift (CrossMap hg19 -> hg38) ---")
with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as fh_in:
    bed_in = fh_in.name
    for _, r in ho.iterrows():
        name = f"{r['gene_symbol']}|{r['circRNA']}"
        fh_in.write(f"{r['hg19_chr']}\t{r['hg19_start']}\t{r['hg19_end']}\t{name}\t0\t{r['strand']}\n")
bed_out = bed_in.replace(".bed", ".hg38.bed")
res = subprocess.run([str(Path.home()/"miniconda3"/"bin"/"CrossMap"), "bed", CHAIN, bed_in, bed_out],
                     capture_output=True, text=True)
print(res.stdout[-800:] if res.stdout else "")
if res.returncode != 0: print("STDERR:", res.stderr[-500:])

lifted_rows = []
with open(bed_out) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6: continue
        chrom, start, end, name, score, strand = p[:6]
        lifted_rows.append({"name": name, "hg38_chr": chrom,
                            "hg38_start": int(start), "hg38_end": int(end),
                            "hg38_strand": strand})
lifted = pd.DataFrame(lifted_rows)
print(f"Lifted {len(lifted)} / {len(ho)} circRNAs ({100*len(lifted)/len(ho):.1f}%)")

# Strip 'chr' prefix to match our universe (Ensembl-style chrom)
lifted["hg38_chr_no"] = lifted["hg38_chr"].str.replace("^chr","", regex=True)
lifted["coord_key"] = (lifted["hg38_chr_no"].astype(str) + ":" +
                       lifted["hg38_start"].astype(str) + "-" +
                       lifted["hg38_end"].astype(str) + ":" + lifted["hg38_strand"])

# Match against our circRNA universe
universe_keys = set(universe["coord_key"])
lifted["in_our_universe"] = lifted["coord_key"].isin(universe_keys)
n_match = int(lifted["in_our_universe"].sum())
print(f"\nExact coordinate matches in our universe: {n_match} / {len(lifted)} = {100*n_match/max(len(lifted),1):.1f}%")

# Match against our BSJ-DE and FSJ-DE sets
bsj_de_keys = set(universe.loc[universe["circRNA"].isin(bsj_add.loc[bsj_add["adj.P.Val"]<0.05, "circRNA"]), "coord_key"])
fsj_de_keys = set(universe.loc[universe["circRNA"].isin(fsj_add.loc[fsj_add["adj.P.Val"]<0.05, "circRNA"]), "coord_key"])
lifted["in_BSJ_DE"] = lifted["coord_key"].isin(bsj_de_keys)
lifted["in_FSJ_DE"] = lifted["coord_key"].isin(fsj_de_keys)

# Hypergeometric: of N=5,820 in our universe, K = Ho-matched (in universe),
# n = our BSJ-DE size (56), k = Ho-matched ∩ BSJ-DE
N = len(universe)
K = n_match
n_bsj = int((bsj_add["adj.P.Val"]<0.05).sum())
k_bsj = int(lifted["in_BSJ_DE"].sum())
P_bsj = hypergeom.sf(k_bsj - 1, N, K, n_bsj)
fold_bsj = (k_bsj/max(n_bsj,1)) / (K/N)
n_fsj = int((fsj_add["adj.P.Val"]<0.05).sum())
k_fsj = int(lifted["in_FSJ_DE"].sum())
P_fsj = hypergeom.sf(k_fsj - 1, N, K, n_fsj)
fold_fsj = (k_fsj/max(n_fsj,1)) / (K/N)
print(f"\n--- Hypergeometric ---")
print(f"  N (our universe):           {N}")
print(f"  K (Ho-sig matched to our universe by coords): {K}")
print(f"  n (our additive BSJ-DE):    {n_bsj}; k (∩ Ho-sig): {k_bsj}; fold {fold_bsj:.2f}; P = {P_bsj:.3g}")
print(f"  n (our additive FSJ-DE):    {n_fsj}; k (∩ Ho-sig): {k_fsj}; fold {fold_fsj:.2f}; P = {P_fsj:.3g}")

# Build a merged table for output
ho_merged = ho.merge(
    lifted.rename(columns={"name":"_name"}).assign(circRNA=lambda d: d["_name"].str.split("|").str[1]),
    on="circRNA", how="left"
).merge(
    universe[["circRNA","coord_key"]].rename(columns={"circRNA":"our_circRNA"}),
    on="coord_key", how="left"
)
ho_merged["our_BSJ_DE_sig"] = ho_merged["our_circRNA"].isin(
    set(bsj_add.loc[bsj_add["adj.P.Val"]<0.05, "circRNA"]))
ho_merged["our_FSJ_DE_sig"] = ho_merged["our_circRNA"].isin(
    set(fsj_add.loc[fsj_add["adj.P.Val"]<0.05, "circRNA"]))
ho_merged.to_csv(OUT_DIR / "Ho2021_circRNAs_vs_ours.tsv", sep="\t", index=False)

# Save the small set of direct hits (in our BSJ-DE)
hits_bsj = ho_merged[ho_merged["our_BSJ_DE_sig"]].copy()
print(f"\n=== Ho circRNAs that are ALSO additive BSJ-DE in our data ===")
print(hits_bsj[["circRNA","gene_symbol","Ho_mean_log2FC","bound","our_circRNA"]].to_string(index=False))
hits_bsj.to_csv(OUT_DIR / "Ho2021_circRNAs_in_our_BSJ_DE.tsv", sep="\t", index=False)

# Summary TSV
pd.DataFrame([{
    "scope": "coord_match_in_universe", "Ho_total": len(ho), "Ho_lifted": len(lifted),
    "Ho_matched_in_universe": K, "our_set_size": N,
    "k": K, "P": np.nan, "fold": np.nan},
    {"scope": "BSJ_DE_intersect", "Ho_total": len(ho), "Ho_lifted": len(lifted),
     "Ho_matched_in_universe": K, "our_set_size": n_bsj,
     "k": k_bsj, "P": P_bsj, "fold": round(fold_bsj, 2)},
    {"scope": "FSJ_DE_intersect", "Ho_total": len(ho), "Ho_lifted": len(lifted),
     "Ho_matched_in_universe": K, "our_set_size": n_fsj,
     "k": k_fsj, "P": P_fsj, "fold": round(fold_fsj, 2)}
]).to_csv(OUT_DIR / "Ho2021_circRNAs_overlap_summary.tsv",
          sep="\t", index=False, float_format="%.4g")
print(f"\nWrote outputs to {OUT_DIR}")
