# Fig S1B — gene PCA (voom log2-CPM), MIDDLE panel of the 3-panel PCA figure.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/05_Summary_Plots/pca_3panel_PSD58.R
# Retains only the gene PCA panel (p2) — PCA of voom log2-CPM from salmon gene
#   counts (filterByExpr + TMM + voom), samples coloured by group, shaped by PSD.
#   The circRNA and CLR PCA panels are dropped.
#
# Strategy (same as Fig_2C / Fig_S1A): keep ALL original plot parameters untouched,
#   render at the original per-panel size, and let ggsave(scale=) squash the figure
#   to 7.75 x 6.83 cm (uniform shrink preserves proportions). Pre-inflate ONLY the
#   axis title/tick fonts AND legend fonts so they land at 8 pt in the output.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(edgeR); library(limma)
  library(ggplot2); library(ggrepel); library(matrixStats)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Supplementary")
GENE_RDS <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_rnaseq_230423/results/star_salmon/salmon.merged.gene_counts.rds"

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

# ---- Gene voom log2-CPM + PCA (identical pipeline to source) ----
se <- readRDS(GENE_RDS); counts <- round(assay(se, "counts"))
samples_X <- paste0("X", samples); mat_g <- counts[, samples_X]; colnames(mat_g) <- samples
group_g <- factor(ifelse(grepl("gHnrnpM", samples), "HnrnpM", "NEG4"),
                  levels = c("NEG4","HnrnpM"))
dge <- DGEList(counts = mat_g, group = group_g)
dge <- dge[filterByExpr(dge, group = group_g), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v   <- voom(dge, model.matrix(~group_g), plot = FALSE)
pca_gene <- do_pca(v$E, sprintf("gene log2-CPM voom (n=%d)", nrow(v$E)))

# ---- Build the gene PCA panel (original parameters preserved) ----
d <- pca_gene
d$group <- factor(ifelse(grepl("gHnrnpM", d$sample), "gHNRNPM", "gCTRL"),
                  levels = c("gCTRL","gHNRNPM"))
d$PSD   <- factor(sub(".*PSD([0-9]+).*", "\\1", d$sample), levels = c("5","8"))
d$short <- gsub("gHnrnpM", "gHNRNPM", gsub("gNEG4", "gCTRL", sub("^[0-9]+_", "", d$sample)))

pc <- ggplot(d, aes(x = PC1, y = PC2, colour = group, shape = PSD, label = short)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(size = 2.25, alpha = 0.9) +
  geom_text_repel(size = 1.3, min.segment.length = 0, seed = 42, show.legend = FALSE,
                  box.padding = 0.3, force = 0.5, force_pull = 5, max.overlaps = Inf) +
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
    axis.title   = element_text(size = axis_render_pt),
    axis.text    = element_text(size = axis_render_pt),
    legend.title = element_text(size = axis_render_pt),
    legend.text  = element_text(size = axis_render_pt),
    plot.margin        = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
    legend.box.spacing = grid::unit(2, "pt"),
    legend.margin      = margin(0, 0, 0, 0)
  )

ggsave(file.path(OUT_DIR, "Fig_S1B.pdf"), pc,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_S1B.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_S1B.pdf)
# ============================================================
# Principal component analysis of gene expression across the seven samples,
# computed on voom log2-CPM values derived from salmon gene counts after filtering
# to expressed genes (edgeR filterByExpr) and TMM normalisation, with zero-variance
# genes removed and features centred and scaled. Each point is one sample,
# positioned by its first two principal components; axis labels report the
# percentage of total variance captured by PC1 and PC2. Points are coloured by
# group, control (gCTRL) and HNRNPM knockdown (gHNRNPM), and shaped by postnatal
# day (circle for day 5, triangle for day 8), with each point labelled by sample.
# Separation of control and knockdown samples indicates a global shift in the gene
# expression programme upon HNRNPM knockdown.
