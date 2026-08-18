# Fig 2A — Per-circRNA mean CLR scatter (gCTRL vs gHNRNPM).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/03_CLR/global_shift_plots_PSD58.R
# Retains only PANEL A of that figure — the per-circRNA mean circular-to-linear
#   ratio (CLR) scatter of control vs HNRNPM knockdown, with the y=x diagonal and
#   BSJ-DE circRNAs highlighted/labelled. Panels B/C/D (density, histogram) dropped.
#
# Strategy: keep ALL original plot parameters (base_size 11, point sizes, repel
#   label size, line widths) untouched, so the plot looks IDENTICAL to the source.
#   Render at the original size and let ggsave(scale=) squash the whole figure down
#   to the target physical size (7.75 x 6.83 cm) — a uniform shrink preserves every
#   proportion. Because that shrink also shrinks the fonts, pre-inflate ONLY the
#   axis title/tick fonts by the same factor so they land at exactly 8 pt.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES     <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
CLR_F   <- file.path(RES, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
BSJ_DE  <- file.path(RES, "01_circRNA_DE", "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")   # additive ~group+PSD (matches Fig 1)
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # final PDF size on the slide
ORIG_W_IN <- 7                     # the source rendered panel A at width = 7 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis title + tick font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Sample groups ----
HM  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
NEG <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")

# ---- Load CLR matrix + BSJ-DE set ----
dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE)
meta_cols <- c("circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id")
mat <- as.matrix(dat[, setdiff(colnames(dat), meta_cols)])
rownames(mat) <- dat$circRNA

bsj_de <- read.table(BSJ_DE, sep = "\t", header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE)
sig_bsj <- bsj_de$circRNA[bsj_de$adj.P.Val < 0.05]

# ---- Per-circRNA mean CLR scatter (ALL original parameters preserved) ----
scat <- data.frame(
  circRNA     = rownames(mat),
  gene_symbol = dat$gene_symbol,
  mean_NEG    = rowMeans(mat[, NEG]),
  mean_HM     = rowMeans(mat[, HM]),
  stringsAsFactors = FALSE
)
scat$sig <- scat$circRNA %in% sig_bsj
scat <- scat[order(scat$sig), ]     # draw sig points on top
n_above   <- sum(scat$mean_HM > scat$mean_NEG)
pct_above <- round(100 * n_above / nrow(scat), 1)

pa <- ggplot(scat, aes(x = mean_NEG, y = mean_HM)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = sig, size = sig, alpha = sig)) +
  geom_text_repel(data = subset(scat, sig), aes(label = gene_symbol),
                  size = 2.7, min.segment.length = 0, max.overlaps = 25, seed = 1) +
  scale_colour_manual(values = c("FALSE" = "grey75", "TRUE" = "#d73027"),
                      labels = c("FALSE" = "NS", "TRUE" = "BSJ-DE (adj.P<0.05)"),
                      name = "") +
  scale_size_manual(values = c("FALSE" = 0.6, "TRUE" = 2.2), guide = "none") +
  scale_alpha_manual(values = c("FALSE" = 0.35, "TRUE" = 0.95), guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Mean CLR per circRNA: gCTRL vs gHNRNPM",
       subtitle = sprintf("%d circRNAs; %d (%.1f%%) above the y=x diagonal",
                          nrow(scat), n_above, pct_above),
       x = "Mean CLR in gCTRL (n=4)",
       y = "Mean CLR in gHNRNPM (n=3)") +
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    legend.position = "bottom",
    axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text  = element_text(size = axis_render_pt)   # -> 8 pt after the ggsave shrink
  )

# scale = SC renders the plot at its original size, then the file is squashed to
# W_CM x H_CM. Uniform shrink => everything keeps the source's proportions.
ggsave(file.path(OUT_DIR, "Fig_2A.pdf"), pa,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_2A.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_2A.pdf)
# ============================================================
# Scatter plot of the mean circular-to-linear ratio (CLR) per circRNA in control
# versus HNRNPM knockdown. Each point is one circRNA. The CLR of a circRNA in a
# sample is 2 times its back-splice junction count divided by the sum of 2 times
# the back-splice junction count and the forward-splice junction count, so it
# ranges from 0 (fully linear) to 1 (fully circular). The horizontal axis shows
# the mean CLR across the four control samples (gCTRL) and the vertical axis shows
# the mean CLR across the three HNRNPM knockdown samples (gHNRNPM), both bounded
# between 0 and 1. The dashed line is the line of identity (y equals x); points
# above it have higher circular fraction in HNRNPM knockdown than in control.
# Red points are circRNAs significantly differentially expressed at the back-splice
# junction level (adjusted P below 0.05 under the additive model that adjusts for
# postnatal day, the same set highlighted in Figure 1) and are labelled with their
# host gene symbols, while grey points are not significant. The subtitle reports
# the total number of circRNAs and the number and percentage lying above the
# diagonal, summarising the transcriptome-wide shift toward higher circular
# fraction under HNRNPM knockdown.
