# HNRNPM regulates circRNA back-splicing in SH-SY5Y neuroblastoma cells

Analysis code for the manuscript *"HNRNPM regulates circRNA back-splicing in SH-SY5Y
neuroblastoma cells"* (submitted to *BMC Genomics*). Voineagu Lab, UNSW Sydney.

Raw and processed RNA-seq data are in GEO under accession
[GSE343178](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343178). HNRNPM eCLIP data
were taken from Ho et al. 2021 (GSE113783).

This repository contains **scripts only**. Every script reads from and writes to absolute paths
on the lab's Katana share, so paths need editing before the code will run elsewhere. The
intermediate result tables are not included; the processed expression tables are published as
Supplementary Tables with the manuscript.

## Layout

```
FinalAnalysis/SCRIPTS/    nine numbered analysis stages, run in order
FinalAnalysis/Figures/    an earlier combined-figure script
CircManuscript/Figures/   per-panel figure scripts (Main and Supplementary)
CircManuscript/d1_*, d2_* circRNA-by-AS overlap and CLR panels
CircManuscript/w0_*..w3_* published-dataset survey and bench-assay statistics
CircManuscript/GEO/       GEO submission packaging
56circ_NeuroCirc/         appends NeuroCirc brain-expression columns to Table S5
```

The analysis stages run in numerical order. `01_circRNA_DE` builds the back-splice junction
(BSJ) and forward-splice junction (FSJ) count matrices from CIRIquant output and runs the
limma-voom differential expression; everything downstream reads its tables.

## Which script made which submitted item

Differential expression uses the additive `~ group + PSD` design throughout, so the tables the
figures read are the `_PSDadd` ones.

| Submitted | Script |
|---|---|
| Fig 1A, 1B | no script; schematic and RT-qPCR chart drawn in PowerPoint |
| Fig 1C | `CircManuscript/Figures/Main/Fig_1A.R` |
| Fig 1D | `CircManuscript/Figures/Main/Fig_1D.R` |
| Fig 2A | `CircManuscript/Figures/Main/Fig_2A.R` |
| Fig 2B | `CircManuscript/Figures/Main/Fig_2Z.R` |
| Fig 2C | `CircManuscript/Figures/Main/Fig_2D.R` |
| Fig 2D | `CircManuscript/Figures/Main/Fig_2C.R` |
| Fig 3A | `CircManuscript/Figures/Main/Fig_3A.R` (the `_v2` output) |
| Fig 3B | `CircManuscript/Figures/Main/Fig_3B.R` |
| Fig 3C | `CircManuscript/Figures/Main/Fig_3C_v2.R` |
| Fig 4A | `CircManuscript/Figures/Main/Fig_1C.R` |
| Fig 4B | `CircManuscript/Figures/Main/Fig_1B_v2.R` |
| Fig 4C | `CircManuscript/Figures/Main/Fig_4F.R` |
| Fig 4D | `CircManuscript/Figures/Supplementary/Fig_S4A_new.R` |
| Fig 5A, 5B | no script; RT-qPCR bars and Sanger chromatograms |
| Fig 6A | `CircManuscript/Figures/Main/Fig_4A.R` |
| Fig 6B | `CircManuscript/Figures/Main/Fig_4E.R` |
| Fig 6C | `CircManuscript/Figures/Main/Fig_4C.R` |
| Fig 6D | `CircManuscript/Figures/Main/Fig_4D.R` |
| Fig 7A, 7B | `CircManuscript/Figures/Main/Fig_5A.R`, `Fig_5B.R` |
| Fig 7C, 7D | `CircManuscript/Figures/Main/Fig_5C.R`, `Fig_5D.R` |
| Fig 7E, 7F | `CircManuscript/Figures/Main/Fig_5E.R`, `Fig_5F.R` |
| Fig S1A, S1B | `CircManuscript/Figures/Supplementary/Fig_S1A.R`, `Fig_S1B.R` |
| Fig S2 | `CircManuscript/Figures/Supplementary/Fig_S3A.R` |
| Fig S3 | `CircManuscript/Figures/Supplementary/Fig_S3C_v1.R` |
| Fig S4 | `CircManuscript/Figures/Supplementary/Fig_S4B_new.R` |
| Fig S5 | `CircManuscript/Figures/Main/Fig_4B.R` |
| Fig S6A-S6D | `CircManuscript/Figures/Supplementary/Fig_S5A.R` to `Fig_S5D.R` |
| Fig S7 | no script; Sanger chromatograms |
| Table S1 | no script; nf-core/rnaseq MultiQC report |
| Table S2 | no script; sample metadata |
| Table S3 | `FinalAnalysis/SCRIPTS/01_circRNA_DE/limma_voom_BSJ_PSD58.R` |
| Table S4 | `FinalAnalysis/SCRIPTS/01_circRNA_DE/FSJ/limma_voom_FSJ_PSD58.R` |
| Table S5 | `FinalAnalysis/SCRIPTS/01_circRNA_DE/FSJ/classify_BSJ_FSJ_additive.py`, then `56circ_NeuroCirc/append_neurocirc_columns.R` |
| Table S6 | `FinalAnalysis/SCRIPTS/02_GeneTx_DE/limma_voom_gene_PSD58.R` |
| Table S7 | `FinalAnalysis/SCRIPTS/02_GeneTx_DE/limma_voom_transcript_PSD58.R` |
| Table S8 | `FinalAnalysis/SCRIPTS/02_GeneTx_DE/enrichment_HnrnpM_PSD58.R` |
| Table S9 | `FinalAnalysis/SCRIPTS/03_CLR/test_global_CLR_shift.py` |
| Table S10 | `FinalAnalysis/SCRIPTS/07_Splicing_rMATS/summarize_rmats_PSD58.R` |
| Table S11 | `CircManuscript/d1_circDE_vs_AS/d1_circDE_vs_AS.R` |
| Table S12 | `FinalAnalysis/SCRIPTS/09_HnrnpM_CLIP_enrichment/test_CLIP_enrichment_PSD58.py` |
| Table S13 | `FinalAnalysis/SCRIPTS/06_Ho2021_Comparison/compare_HnrnpM_DE_to_Ho2021_PSD58.py` |
| Table S14, S15 | no script; primer and gRNA sequences |

### The script names predate the final figure numbering

The scripts were written against an earlier five-figure layout and were never renamed, so a
script's name does not always say where its panel ended up. `Fig_1B_v2.R` and `Fig_1C.R` became
Figure 4B and 4A. `Fig_4A/4C/4D/4E.R` became Figure 6A/6C/6D/6B, and `Fig_4F.R` became Figure 4C.
`Fig_5A`-`Fig_5F.R` became Figure 7A-7F. `Fig_4B.R` sits in `Figures/Main` but was submitted as
Supplementary Figure S5. Within Figure 2 the panel letters also swap, and Supplementary
`Fig_S3A.R` and `Fig_S3C_v1.R` became Supplementary Figures S2 and S3.

### Scripts that did not feed the manuscript

Also included, for completeness, are the FSJ threshold-selection comparison, the
`05_Summary_Plots` and other earlier combined-panel scripts, superseded figure variants
(`Fig_1B`, `Fig_1E`, `Fig_2B`, `Fig_3C`, `Fig_3C_v1`, `Fig_3A_nonadditive`, `Fig_S3C`,
`Fig_S4A`, `Fig_S4B`, `Fig_Sz`), `d1_circDE_vs_AS_any.R`, `d2_CLR_violins.R`, the GEO packaging
scripts, and the `w0_circ_Disease` survey of published circRNA disease datasets. Anything not
listed in the table above did not produce submitted material.

## Software

Read alignment and quantification used nf-core/rnaseq 3.11.1 and nf-core/circrna under
Nextflow 22.10.7 on the UNSW Katana HPC, with STAR 2.7.9a.
circRNAs were called as the intersection of CIRCexplorer2 and CIRIquant 1.1.2, and BSJ and FSJ
counts were taken from the CIRIquant GTF attributes. Differential alternative splicing used
rMATS-turbo 4.3.0 (JCEC). Downstream analysis used R 4.3.3 with limma, edgeR, clusterProfiler
and ggplot2 3.5.1, and Python 3 with pandas, numpy, scipy and pyBigWig.
