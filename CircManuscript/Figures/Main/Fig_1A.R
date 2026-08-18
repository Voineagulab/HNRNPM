# Fig 1A — BSJ volcano (labelled) for the HnrnpM circRNA manuscript.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/01_circRNA_DE/limma_voom_BSJ_PSD58.R
# Draws the BSJ DE volcano (additive model, ~ group + PSD) and labels the 10
#   most significant DE circRNAs (top 10 by adj.P.Val; adapted from
#   ForPublication/FinalAnalysis/Figures/Figure1_circRNA_DE.R).
# Saved as PDF sized for the manuscript slide (W_CM x H_CM), fonts >= 8 pt where
#   set, via cairo_pdf, to Figures/Main/Fig_1A.pdf.
# CIRIquant BSJ counts, 7 samples (3 HnrnpM_PSD5/8 vs 4 NEG4_PSD5/8),
# library-size normalisation from MultiQC R1 reads.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(limma); library(edgeR); library(ggplot2); library(ggrepel)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
COUNTS_F <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/bsj_matrix_CIRIquant_HnrnpM_PSD58.tsv")
MULTIQC  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_circrna_230420/results/multiqc/multiqc_data/multiqc_general_stats.txt"
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

ANNOT_COLS <- c("gene_symbol","chr","start","end","strand","circ_type","gene_id")

# ---- Target physical size on the slide ----
W_CM   <- 5.24; H_CM <- 6.83      # width x height as it should appear in PowerPoint
MIN_PT   <- 8               # base font size for axes (points)
LABEL_PT <- 4               # gene-label font size (points) — lower = smaller labels
TITLE_PT <- 3               # title + subtitle font size (points)
LEG_PT   <- 3               # legend text + title font size (points)
lab_size <- LABEL_PT / .pt  # geom_text_repel size (mm) for LABEL_PT

# ---- Load counts + annotation ----
raw <- read.table(COUNTS_F, header = TRUE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE)
annot   <- raw[, c("circRNA", ANNOT_COLS), drop = FALSE]
counts  <- round(as.matrix(raw[, !(colnames(raw) %in% c("circRNA", ANNOT_COLS))]))
rownames(counts) <- raw$circRNA
samples <- colnames(counts)
cat("Loaded counts:", nrow(counts), "circRNAs x", ncol(counts), "samples\n")

group <- factor(ifelse(grepl("gHnrnpM", samples), "gHnrnpM", "gNEG4"),
                levels = c("gNEG4","gHnrnpM"))
PSD   <- factor(sub(".*PSD([0-9]+).*", "\\1", samples), levels = c("5","8"))

# ---- Library-size from MultiQC R1 ----
stats <- read.table(MULTIQC, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
r1    <- stats[grepl("_1$", rownames(stats)),
               "FastQC_mqc-generalstats-fastqc-total_sequences", drop = FALSE]
rownames(r1) <- sub("_1$", "", rownames(r1))
lib_sizes <- r1[samples, 1]; names(lib_sizes) <- samples

dge <- DGEList(counts = counts, lib.size = lib_sizes, group = group)

# ---- Additive model: ~ group + PSD (this is the panel model) ----
design_add <- model.matrix(~ group + PSD); rownames(design_add) <- samples
v   <- voom(dge, design_add, plot = FALSE)
fit <- eBayes(lmFit(v, design_add))
res <- topTable(fit, coef = "groupgHnrnpM", number = Inf, sort.by = "P")
res$gene_symbol <- annot$gene_symbol[match(rownames(res), annot$circRNA)]

# ---- Labelled volcano (reproduces panel_BSJ_volcano.png) ----
res$direction <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0, "UP in HNRNPM KD",
                 ifelse(res$adj.P.Val < 0.05 & res$logFC < 0, "DOWN in HNRNPM KD", "NS"))
n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
n_up  <- sum(res$direction == "UP in HNRNPM KD")
n_dn  <- sum(res$direction == "DOWN in HNRNPM KD")
sig_pts <- subset(res, direction != "NS")               # all DE circRNAs (adj.P<0.05)
# Label only the 10 most significant DE circRNAs (adapted from
# ForPublication/FinalAnalysis/Figures/Figure1_circRNA_DE.R lines 75-76):
sig_pts <- sig_pts[order(sig_pts$adj.P.Val), ][1:10, ]  # top 10 by adj.P.Val
cat(sprintf("tested %d, sig %d (UP %d / DOWN %d); labelling top %d by adj.P.Val\n",
            nrow(res), n_sig, n_up, n_dn, nrow(sig_pts)))

pv <- ggplot(res, aes(x = logFC, y = -log10(P.Value), colour = direction)) +
  geom_point(alpha = 0.5, size = 0.15) +                                         # JW: dot size
  scale_colour_manual(values = c("UP in HNRNPM KD"="#d73027",
                                 "DOWN in HNRNPM KD"="#2166ac","NS"="grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  labs(title = "limma-voom BSJ: ~group + PSD  (FinalAnalysis 01)",
       subtitle = sprintf("%d tested  ·  %d sig (UP %d / DOWN %d) at adj.P<0.05",
                          nrow(res), n_sig, n_up, n_dn),
       x = "log2 Fold Change", y = "-log10(P)", colour = "") +
  theme_bw(base_size = MIN_PT) +
  theme(
    legend.position = "bottom",
    text          = element_text(size = MIN_PT),
    axis.text     = element_text(size = MIN_PT),
    axis.title    = element_text(size = MIN_PT),
    # --- JW: title/subtitle take less vertical space ---
    plot.title    = element_text(size = TITLE_PT, hjust = 0.5, margin = margin(b = 1)),
    plot.subtitle = element_text(size = TITLE_PT, hjust = 0.5, margin = margin(b = 1)),
    # --- JW: compact legend ---
    legend.text        = element_text(size = LEG_PT),
    legend.title       = element_blank(),              # empty title — reclaim its space
    legend.key.size    = grid::unit(6, "pt"),          # smaller colour keys
    legend.margin      = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(1, "pt"),          # gap between plot and legend
    legend.spacing     = grid::unit(2, "pt")
  ) +
  ggrepel::geom_text_repel(
    data = sig_pts, aes(label = gene_symbol, colour = direction),
    size = lab_size,
    fontface = "italic",
    box.padding = 0.05,       # JW tighter: minimal margin around each label
    point.padding = 0.12,     # JW small clear gap around the dot (keeps label off it, but close)
    force = 0.3,              # JW low mutual repulsion between labels
    force_pull = 2,           # JW stronger pull: labels hug their dots
    min.segment.length = Inf,   # JW: never draw a connector line
    segment.color = NA,         # JW: no line even if a label is nudged
    max.iter = 100000,
    max.overlaps = Inf,
    seed = 1,
    show.legend = FALSE)

ggsave(file.path(OUT_DIR, "Fig_1A.pdf"), pv,
       width = W_CM, height = H_CM, units = "cm",
       device = cairo_pdf)


# ---- Figure legend ----
# Figure 1. Identification of HNRNPM-dependent circRNAs in SH-SY5Y cells.
# (A) Volcano plot of the back-splice junction (BSJ) fold change (FC) in 5820 of circRNAs detected in SH-SY5Y cells treated with HNRNPM-targeting gRNA (gHNRNPM) or control gRNAs (gCTRL). The 55 red and 1 blue dots are the differentially expressed (DE) circRNA with significantly higher and lower, respectively, expression in the HNRNPM knockdown (KD) samples as tested with limma-voom on CIRIquant BSJ read counts under an additive design that adjusts for cells harvested at different days post gRNA-selection (PSD) across three KD and four control samples. Only circRNAs with an adjusted P value < 0.05 are coloured. The horizontal and vertical dashed line marks the  P value < 0.05 and |log2(FC)| > 1 threshold, respectively. 
