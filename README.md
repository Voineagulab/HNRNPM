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

## Directory structure

The analysis stages run in numerical order, and the scripts assume the structure below. Stage
`01_circRNA_DE` builds the back-splice junction (BSJ) and forward-splice junction (FSJ) count
matrices from CIRIquant output and runs the differential expression; everything downstream reads
its tables. Differential expression uses the additive `~ group + PSD` design throughout, so the
tables the figure scripts read are the ones suffixed `_PSDadd`.

```
HNRNPM
├── FinalAnalysis
│   ├── SCRIPTS
│   │   ├── 01_circRNA_DE               circRNA quantification, BSJ differential expression
│   │   │   └── FSJ                     cognate linear junction DE, BSJ x FSJ classification
│   │   │       └── ThresholdSelection  sensitivity of the classification to the FSJ threshold
│   │   ├── 02_GeneTx_DE                gene and transcript DE, GO over-representation
│   │   ├── 03_CLR                      circular-to-linear ratio and its global shift
│   │   ├── 05_Summary_Plots            combined classification and PCA overview panels
│   │   ├── 06_Ho2021_Comparison        replication against HNRNPM knockdown in PC-3M cells
│   │   ├── 07_Splicing_rMATS           differential alternative splicing
│   │   ├── 08_circRNA_AS_overlap       circRNA and AS event splice-site co-localisation
│   │   └── 09_HnrnpM_CLIP_enrichment   HNRNPM eCLIP signal at circRNA-flanking introns
│   └── Figures                         combined figure script for the circRNA DE analysis
├── CircManuscript
│   ├── Figures
│   │   ├── Main                        per-panel scripts for the main figures
│   │   └── Supplementary               per-panel scripts for the supplementary figures
│   ├── d1_circDE_vs_AS                 circRNA x differential-AS overlap, Fisher's exact tests
│   ├── d2_CLR                          circular-to-linear ratio distribution panels
│   ├── w0_circ_Disease                 survey of published circRNA disease datasets
│   │   └── SCRIPTS
│   ├── w1_EdU                          EdU S-phase fraction statistics
│   ├── w2_growth_curve                 growth curves, ANOVA and doubling times
│   ├── w3_circ_FC_ceiling              proliferation-only circRNA fold-change ceiling
│   └── GEO                             packaging of the GEO submission
│       ├── md5_checksums
│       └── processed_files
└── 56circ_NeuroCirc                    NeuroCirc brain-expression annotation
```

## Software

Read alignment and quantification used nf-core/rnaseq 3.11.1 and nf-core/circrna under
Nextflow 22.10.7 on the UNSW Katana HPC, with STAR 2.7.9a.
circRNAs were called as the intersection of CIRCexplorer2 and CIRIquant 1.1.2, and BSJ and FSJ
counts were taken from the CIRIquant GTF attributes. Differential alternative splicing used
rMATS-turbo 4.3.0 (JCEC). Downstream analysis used R 4.3.3 with limma, edgeR, clusterProfiler
and ggplot2 3.5.1, and Python 3 with pandas, numpy, scipy and pyBigWig.
