"""
FinalAnalysis item 09 step 1 — re-derive flanking-intron annotation against
the new CIRIquant universe (5,820 circRNAs) and the additive gene-level DE
(publication primary).

Outputs:
  flanking_introns_per_circRNA.tsv     (per-circRNA 5'/3' flank coords)
  within_gene_background_regions.tsv   (~5 per circRNA, length-matched)
  non_host_gene_background_regions.tsv (5,000 regions from expressed non-host genes)
"""
import re
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV")
GTF  = "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/DATA/Homo_sapiens.GRCh38.109.gtf"
BSJ  = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "01_circRNA_DE" / "bsj_matrix_CIRIquant_HnrnpM_PSD58.tsv"
LIMMA_GENE_ADD = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "02_GeneTx_DE" / "limma_gene_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"
OUT_DIR = ROOT / "ForPublication" / "FinalAnalysis" / "RESULTS" / "09_HnrnpM_CLIP_enrichment"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_BG_PER_CIRC = 5
N_BG_NON_HOST = 5000

# ---- Load circRNAs from the new 5,820 universe ----
mat = pd.read_csv(BSJ, sep="\t")
mat["chr"]    = mat["chr"].astype(str)
mat["start"]  = mat["start"].astype(int)
mat["end"]    = mat["end"].astype(int)
mat["strand"] = mat["strand"].astype(str)
mat["primary_gene"] = mat["gene_id"].apply(lambda s: str(s).split(",")[0].strip())
print(f"circRNAs in universe: {len(mat)}")

counts = mat.set_index("circRNA")
host_genes = set(counts["primary_gene"])
print(f"Unique host genes:    {len(host_genes)}")

# Expressed-gene universe = limma additive-tested set (gene level)
lm = pd.read_csv(LIMMA_GENE_ADD, sep="\t")
expressed_genes = set(lm["gene_id"].astype(str))
print(f"Expressed genes (additive limma universe): {len(expressed_genes)}")

# ---- Parse GTF: collect exons per gene ----
print("\nParsing GTF (~30s) ...")
gene_id_re = re.compile(r'gene_id "([^"]+)"')
exons_by_gene = {}
n_exon_lines = 0
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "exon": continue
        m = gene_id_re.search(f[8])
        if not m: continue
        g = m.group(1)
        exons_by_gene.setdefault(g, []).append((f[0], int(f[3])-1, int(f[4]), f[6]))
        n_exon_lines += 1
print(f"  exon GTF lines: {n_exon_lines}; unique genes: {len(exons_by_gene)}")

def merge_intervals(ivals):
    ivals = sorted(ivals)
    merged = []
    cur_s, cur_e = ivals[0]
    for s, e in ivals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e)); cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged

introns_by_gene = {}
for g, exons in exons_by_gene.items():
    chrom = exons[0][0]; strand = exons[0][3]
    se = [(s, e) for (_, s, e, _) in exons]
    merged = merge_intervals(se)
    if len(merged) < 2: continue
    intr = []
    for i in range(len(merged) - 1):
        s0, e0 = merged[i][1], merged[i+1][0]
        if e0 > s0:
            intr.append((chrom, s0, e0, strand))
    if intr:
        introns_by_gene[g] = intr
print(f"  genes with >=1 intron: {len(introns_by_gene)}")

# ---- Per-circRNA 5' and 3' flanks ----
flank_rows = []
for circ_id, r in counts.iterrows():
    g = r["primary_gene"]
    introns = introns_by_gene.get(g, [])
    flank5 = None; flank3 = None
    for (chrom, s, e, strand) in introns:
        if chrom != r["chr"]: continue
        if e == r["start"]: flank5 = (chrom, s, e, strand)
        if s == r["end"]:   flank3 = (chrom, s, e, strand)
    flank_rows.append({
        "circRNA":     circ_id,
        "gene_symbol": r["gene_symbol"],
        "gene_id":     g,
        "chr":         r["chr"],
        "strand":      r["strand"],
        "circ_start":  r["start"],
        "circ_end":    r["end"],
        "flank5_start":  flank5[1] if flank5 else np.nan,
        "flank5_end":    flank5[2] if flank5 else np.nan,
        "flank5_length": (flank5[2]-flank5[1]) if flank5 else np.nan,
        "flank3_start":  flank3[1] if flank3 else np.nan,
        "flank3_end":    flank3[2] if flank3 else np.nan,
        "flank3_length": (flank3[2]-flank3[1]) if flank3 else np.nan,
        "n_flanks_found": int(flank5 is not None) + int(flank3 is not None)
    })
flank_df = pd.DataFrame(flank_rows)
flank_df.to_csv(OUT_DIR / "flanking_introns_per_circRNA.tsv", sep="\t", index=False)
print(f"\nFlanking-intron resolution:")
print(flank_df["n_flanks_found"].value_counts().to_string())
print(f"  >=1 flank found: {(flank_df['n_flanks_found']>=1).sum()} / {len(flank_df)} = {100*(flank_df['n_flanks_found']>=1).mean():.1f}%")
print(f"  Both flanks:    {(flank_df['n_flanks_found']==2).sum()}")

# ---- Within-host-gene background regions ----
within_rows = []
for _, r in flank_df.iterrows():
    if r["n_flanks_found"] == 0: continue
    g = r["gene_id"]
    introns = introns_by_gene.get(g, [])
    if len(introns) < 2: continue
    target_len = int(np.nanmean([r["flank5_length"], r["flank3_length"]]))
    if target_len < 50: target_len = 50
    non_flank = []
    for (chrom, s, e, strand) in introns:
        if chrom != r["chr"]: continue
        if (s == r["flank5_start"] and e == r["flank5_end"]): continue
        if (s == r["flank3_start"] and e == r["flank3_end"]): continue
        non_flank.append((chrom, s, e, strand))
    if not non_flank: continue
    for _ in range(N_BG_PER_CIRC):
        i = RNG.integers(0, len(non_flank))
        ch, s, e, st = non_flank[i]
        L = min(target_len, e - s)
        if L < 50: continue
        if e - s == L:
            rs, re_ = s, e
        else:
            rs = int(RNG.integers(s, e - L + 1)); re_ = rs + L
        within_rows.append({"circRNA": r["circRNA"], "gene_id": g,
                            "chr": ch, "start": int(rs), "end": int(re_), "strand": st,
                            "length": int(re_ - rs)})
within_df = pd.DataFrame(within_rows)
within_df.to_csv(OUT_DIR / "within_gene_background_regions.tsv", sep="\t", index=False)
print(f"\nWithin-gene background regions: {len(within_df)} (across {within_df['circRNA'].nunique()} circRNAs)")

# ---- Non-host-gene background pool ----
non_host_genes = expressed_genes - host_genes
target_median_len = int(np.nanmedian(np.concatenate([
    flank_df["flank5_length"].dropna().values,
    flank_df["flank3_length"].dropna().values
])))
print(f"\nMedian flanking-intron length: {target_median_len:.0f} bp")
print(f"Non-host expressed genes pool: {len(non_host_genes)}")

non_host_intr = []
for g in non_host_genes:
    for (chrom, s, e, strand) in introns_by_gene.get(g, []):
        if e - s >= target_median_len:
            non_host_intr.append((g, chrom, s, e, strand))
print(f"Non-host introns >= {target_median_len} bp: {len(non_host_intr)}")

non_host_rows = []
indices = RNG.integers(0, len(non_host_intr), size=N_BG_NON_HOST)
for i in indices:
    g, ch, s, e, st = non_host_intr[i]
    L = target_median_len
    if e - s == L:
        rs, re_ = s, e
    else:
        rs = int(RNG.integers(s, e - L + 1)); re_ = rs + L
    non_host_rows.append({"gene_id": g, "chr": ch, "start": int(rs),
                          "end": int(re_), "strand": st, "length": int(re_ - rs)})
non_host_df = pd.DataFrame(non_host_rows)
non_host_df.to_csv(OUT_DIR / "non_host_gene_background_regions.tsv", sep="\t", index=False)
print(f"Non-host-gene background regions sampled: {len(non_host_df)}")
print("\nDone.")
