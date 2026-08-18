# Fig 2C — CLR PCA (gCTRL vs gHNRNPM), rightmost panel of the 3-panel PCA figure.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/05_Summary_Plots/pca_3panel_PSD58.R
# Retains only the CLR PCA panel (p3) of that figure — PCA of the per-circRNA,
#   per-sample CLR matrix (5,820 circRNAs), samples coloured by group and shaped by
#   PSD. The circRNA-CPM and gene-voom PCA panels are dropped.
#
# Strategy (same as Fig_2A/2B/2D): keep ALL original plot parameters (base_size 11,
#   point size, repel label size, guide lines) untouched so the plot looks IDENTICAL
#   to the source panel. Render at the original per-panel size and let ggsave(scale=)
#   squash the whole figure down to 7.75 x 6.83 cm — a uniform shrink preserves every
#   proportion. Pre-inflate ONLY the axis title/tick fonts AND legend fonts by the
#   same factor so they land at exactly 8 pt in the output.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(matrixStats)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES     <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
CLR_F   <- file.path(RES, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # final PDF size on the slide
ORIG_W_IN <- 13/3                  # each panel of the source 3-panel figure (13 in wide) is ~13/3 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Sample order (controls first, then knockdown) ----
HM  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
NEG <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
samples <- c(NEG, HM)

# ---- PCA helper (identical to source) ----
do_pca <- function(mat, label) {
  rv <- rowVars(mat); mat <- mat[!is.na(rv) & rv > 0, , drop = FALSE]
  pr <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  v <- 100 * pr$sdev^2 / sum(pr$sdev^2)
  data.frame(sample = colnames(mat), PC1 = pr$x[,1], PC2 = pr$x[,2],
             label = label, pc1 = round(v[1],1), pc2 = round(v[2],1),
             n_features = nrow(mat))
}

# ---- Load CLR matrix + PCA ----
clr_dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
clr_mat <- as.matrix(clr_dat[, samples]); rownames(clr_mat) <- clr_dat$circRNA
pca_clr <- do_pca(clr_mat, sprintf("CLR (%d circRNAs)", nrow(clr_mat)))

# ---- Build the CLR PCA panel (original parameters preserved) ----
d <- pca_clr
d$group <- factor(ifelse(grepl("gHnrnpM", d$sample), "gHNRNPM", "gCTRL"),
                  levels = c("gCTRL","gHNRNPM"))
d$PSD   <- factor(sub(".*PSD([0-9]+).*", "\\1", d$sample), levels = c("5","8"))
# Point labels: strip the leading "N_" and relabel groups gNEG4 -> gCTRL, gHnrnpM -> gHNRNPM
d$short <- gsub("gHnrnpM", "gHNRNPM", gsub("gNEG4", "gCTRL", sub("^[0-9]+_", "", d$sample)))

pc <- ggplot(d, aes(x = PC1, y = PC2, colour = group, shape = PSD, label = short)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(size = 2.25, alpha = 0.9) +
  geom_text_repel(size = 1.3, min.segment.length = 0, seed = 42, show.legend = FALSE,
                  box.padding = 0.3, force = 0.5, force_pull = 5, max.overlaps = Inf) +  # strong pull to points, avoid drift
  scale_colour_manual(values = c(gCTRL = "#4C78A8", gHNRNPM = "#E45756")) +
  scale_shape_manual(values = c(`5` = 16, `8` = 17)) +
  labs(title = unique(d$label),
       subtitle = sprintf("%d features", unique(d$n_features)),
       x = sprintf("PC1 (%.1f%%)", unique(d$pc1)),
       y = sprintf("PC2 (%.1f%%)", unique(d$pc2)),
       colour = "Group", shape = "PSD") +
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    legend.position = "right",
    axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text  = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    plot.margin        = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),  # trim outer margins (esp. right)
    legend.box.spacing = grid::unit(2, "pt"),                              # shrink gap between plot and legend
    legend.margin      = margin(0, 0, 0, 0)
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_2C.pdf"), pc,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_2C.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_2C.pdf)
# ============================================================
# Principal component analysis of the per-circRNA circular-to-linear ratio (CLR)
# across the seven samples. The CLR of a circRNA in a sample is 2 times its
# back-splice junction count divided by the sum of 2 times the back-splice junction
# count and the forward-splice junction count. PCA was computed on the CLR matrix
# of all detected circRNAs after removing zero-variance rows, with features centred
# and scaled. Each point is one sample, positioned by its first two principal
# components; the axis labels report the percentage of total variance captured by
# PC1 and PC2. Points are coloured by group, control (gCTRL) and HNRNPM knockdown
# (gHNRNPM), and shaped by postnatal day (circle for day 5, triangle for day 8),
# with each point labelled by sample. Separation of the control and knockdown
# samples along the principal components indicates that HNRNPM knockdown produces a
# coordinated, sample-wide change in circular-to-linear stoichiometry.
