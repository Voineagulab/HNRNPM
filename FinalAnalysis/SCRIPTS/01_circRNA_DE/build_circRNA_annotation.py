"""
Generate circRNA_annotation.tsv (reproducibility script).

For every circRNA detected by CIRCexplorer2 in any of the 24 nf-core/circrna
samples, record:
  - circRNA   : BED-style ID (chr:start-end:strand)
  - gene_symbol : gene name of the first host gene (Ensembl 109 GTF lookup;
                  NA for intergenic circRNAs); single string (no commas)
  - gene_id   : Ensembl gene_id(s) from CIRCexplorer2 annotation, comma-separated
                for multi-host circRNAs (e.g. nested or overlapping genes)
  - circ_type : circRNA / ciRNA / EIciRNA, from CIRCexplorer2 annotation

CIRCexplorer2 BED layout follows the nf-core/circrna `nfcore_Output_columnNames`
spec; columns 0..15 are documented in `DATA/nfcore_Output_columnNames.csv`.

Multi-host gene_symbol convention: take the FIRST ENSG in the (possibly
comma-separated) gene_id field; matches the bucket-classification convention
in CLAUDE.md (Session 11 step 31).
"""
import os
import re
from pathlib import Path

BASE = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_circrna_230420/results/circrna_discovery/circexplorer2")
GTF  = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/DATA/Homo_sapiens.GRCh38.109.gtf")
OUT  = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/circRNA_annotation.tsv")
OUT.parent.mkdir(parents=True, exist_ok=True)

# CE2 BED column indices (0-based) — matches DATA/nfcore_Output_columnNames.csv
COL_CIRC_ID  = 3   # "name"
COL_CIRC_TYPE = 12  # "circType"
COL_GENE_ID  = 13  # "geneEnsemblName" (parent gene ID(s), comma-separated)

# ---- Step 1: parse Ensembl 109 GTF for ENSG -> gene_name map ----
print(f"Parsing GTF: {GTF}")
ensg_to_symbol = {}
re_gene = re.compile(r'gene_id "([^"]+)"')
re_name = re.compile(r'gene_name "([^"]+)"')
n_genes = 0
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] != "gene":
            continue
        attrs = cols[8]
        m_gene = re_gene.search(attrs)
        m_name = re_name.search(attrs)
        if not m_gene:
            continue
        ensg = m_gene.group(1)
        symbol = m_name.group(1) if m_name else "NA"
        ensg_to_symbol[ensg] = symbol
        n_genes += 1
print(f"  Parsed {n_genes} gene records ({len(ensg_to_symbol)} unique ENSG -> symbol)")

# ---- Step 2: collect (circRNA, gene_id, circ_type) from all CE2 BEDs ----
samples = sorted([d.name for d in BASE.iterdir() if d.is_dir() and d.name != "intermediates"])
print(f"\nFound {len(samples)} CIRCexplorer2 samples to parse")

records = {}   # circRNA -> (gene_id_str, circ_type)
n_seen = 0
n_unique = 0
for s in samples:
    bed = BASE / s / f"{s}.bed"
    if not bed.exists():
        print(f"  [warn] missing: {bed}")
        continue
    with open(bed) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= COL_GENE_ID:
                continue
            cid = p[COL_CIRC_ID]
            gid = p[COL_GENE_ID] if p[COL_GENE_ID] else "NA"
            ctype = p[COL_CIRC_TYPE] if len(p) > COL_CIRC_TYPE else "NA"
            n_seen += 1
            if cid not in records:
                records[cid] = (gid, ctype)
                n_unique += 1
print(f"  Total CE2 BED rows read: {n_seen:,}")
print(f"  Unique circRNAs across 24 samples: {n_unique:,}")

# ---- Step 3: map first ENSG in gene_id -> gene_symbol ----
def first_symbol(gid):
    if not gid or gid == "NA":
        return "NA"
    first = gid.split(",")[0]
    return ensg_to_symbol.get(first, "NA")

# ---- Step 4: write TSV ----
n_intergenic = 0
n_multi_host = 0
with open(OUT, "w") as f:
    f.write("circRNA\tgene_symbol\tgene_id\tcirc_type\n")
    for cid in sorted(records.keys()):
        gid, ctype = records[cid]
        if gid == "NA" or gid == "":
            sym = "NA"
            n_intergenic += 1
        else:
            sym = first_symbol(gid)
            if "," in gid:
                n_multi_host += 1
            if sym == "NA":
                n_intergenic += 1   # ENSG present but unresolved against this GTF
        f.write(f"{cid}\t{sym}\t{gid}\t{ctype}\n")

print(f"\nWrote: {OUT}")
print(f"  Total rows:           {len(records):,}")
print(f"  Intergenic / NA symbol: {n_intergenic:,} ({100*n_intergenic/len(records):.1f}%)")
print(f"  Multi-host gene_id (comma in gene_id): {n_multi_host:,} ({100*n_multi_host/len(records):.1f}%)")

# ---- Compare to existing table (sanity check) ----
EXISTING = Path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/Preliminary/circRNA_annotation.tsv")
if EXISTING.exists():
    import pandas as pd
    old = pd.read_csv(EXISTING, sep="\t")
    new = pd.read_csv(OUT, sep="\t")
    print(f"\n--- Comparison to existing {EXISTING.name} ---")
    print(f"  Existing rows: {len(old):,}")
    print(f"  New rows:      {len(new):,}")
    shared = set(old["circRNA"]) & set(new["circRNA"])
    old_only = set(old["circRNA"]) - shared
    new_only = set(new["circRNA"]) - shared
    print(f"  Shared circRNA IDs: {len(shared):,}")
    print(f"  Existing-only:      {len(old_only):,}")
    print(f"  New-only:           {len(new_only):,}")
    # For shared, how often do gene_symbol / gene_id / circ_type agree?
    old_idx = old.set_index("circRNA")
    new_idx = new.set_index("circRNA")
    shared_l = sorted(shared)
    agree_sym = sum(str(old_idx.loc[c, "gene_symbol"]) == str(new_idx.loc[c, "gene_symbol"]) for c in shared_l)
    agree_gid = sum(str(old_idx.loc[c, "gene_id"])     == str(new_idx.loc[c, "gene_id"])     for c in shared_l)
    agree_typ = sum(str(old_idx.loc[c, "circ_type"])   == str(new_idx.loc[c, "circ_type"])   for c in shared_l)
    print(f"  Agreement on shared IDs:")
    print(f"    gene_symbol: {agree_sym:,} / {len(shared_l):,} ({100*agree_sym/len(shared_l):.1f}%)")
    print(f"    gene_id:     {agree_gid:,} / {len(shared_l):,} ({100*agree_gid/len(shared_l):.1f}%)")
    print(f"    circ_type:   {agree_typ:,} / {len(shared_l):,} ({100*agree_typ/len(shared_l):.1f}%)")
