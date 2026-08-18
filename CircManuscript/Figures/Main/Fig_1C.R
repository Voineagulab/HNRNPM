# Fig 1C — HnrnpM eCLIP enrichment at circRNA flanking introns (violin plot only).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/09_HnrnpM_CLIP_enrichment/plot_enrichment_PSD58.R
# Retains only PANEL B of that figure — the per-circRNA mean flank-score violin
#   plot (the bottom panel; the upper per-region box panel A is dropped).
# DE = additive BSJ-DE circRNAs (~ group + PSD, adj.P<0.05, publication primary).
# Saved as PDF sized for the manuscript slide (W_CM x H_CM), fonts 8 pt, via
#   cairo_pdf, to Figures/Main/Fig_1C.pdf.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES09   <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/09_HnrnpM_CLIP_enrichment")
SCORE   <- file.path(RES09, "circRNA_flank_score_PSD58.tsv")
TESTS   <- file.path(RES09, "enrichment_test_summary.tsv")
BSJ_ADD <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size on the slide ----
W_CM   <- 4.09; H_CM <- 6.83      # width x height as it should appear in PowerPoint
MIN_PT <- 8                        # font size (points)

# ---- Load data ----
sc    <- read.table(SCORE,   sep = "\t", header = TRUE, check.names = FALSE)
tests <- read.table(TESTS,   sep = "\t", header = TRUE, check.names = FALSE)
bsj   <- read.table(BSJ_ADD, sep = "\t", header = TRUE, check.names = FALSE)

# Additive BSJ-DE circRNAs (publication primary)
de_circs <- bsj$circRNA[bsj$adj.P.Val < 0.05]
de_n  <- sum(sc$circRNA %in% de_circs)
nde_n <- sum(!(sc$circRNA %in% de_circs))

# ---- Per-circRNA mean flank score violin (panel B of the source figure) ----
sc$is_DE_lab <- ifelse(sc$circRNA %in% de_circs,
                       sprintf("DE \n(n=%d)", de_n),
                       sprintf("non-DE \n(n=%d)", nde_n))
p_de_lab  <- sprintf("DE \n(n=%d)", de_n)
p_nde_lab <- sprintf("non-DE \n(n=%d)", nde_n)
p_anno <- tests$p_one_sided[tests$DE_axis == "BSJ_additive" &
                              tests$test == "1. DE flank vs non-DE flank"]

pc <- ggplot(sc, aes(x = is_DE_lab, y = mean_flank_log2, fill = is_DE_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_violin(alpha = 0.5, colour = "grey20", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  scale_fill_manual(values = setNames(c("#d73027","#4C78A8"), c(p_de_lab, p_nde_lab)),
                    guide = "none") +
  labs(title = "Per-circRNA mean flank score (DE vs non-DE)",
       subtitle = sprintf("Wilcoxon rank-sum (one-sided): P = %.2e", p_anno),
       x = NULL, y = "Mean log2(IP/Input) at flank introns") +
  theme_bw(base_size = MIN_PT) +
  theme(
    text          = element_text(size = MIN_PT),
    axis.text     = element_text(size = MIN_PT),
    axis.title    = element_text(size = MIN_PT),
    plot.title    = element_text(size = 3, hjust = 0.5, margin = margin(b = 1)),
    plot.subtitle = element_text(size = 4, hjust = 0, margin = margin(b = 1))
  )

ggsave(file.path(OUT_DIR, "Fig_1C.pdf"), pc,
       width = W_CM, height = H_CM, units = "cm",
       device = cairo_pdf)
cat("Wrote:", file.path(OUT_DIR, "Fig_1C.pdf"), "\n")

# # ---- Figure legend ----
# (C) Violin and box plot comparing the per-circRNA mean HNRNPM eCLIP flank score between DE and non-DE circRNAs. For each circRNA the flank score is the mean log2 of IP/Input across its flanking introns. The horizontal axis separates DE circRNAs from non-DE circRNAs, with the number of circRNAs in each group shown in parentheses. The vertical axis shows the mean log2 of IP/Input at flanking introns. Violins show the full distribution, the inset white box shows the median and interquartile range, and the horizontal dashed line marks a log2 value of zero. The annotation reports the one-sided Wilcoxon rank sum P value testing whether DE circRNAs have higher flank binding than non-DE circRNAs.