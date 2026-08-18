# Fig 4A — Gene-level volcano (additive ~ group + PSD, publication primary).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/02_GeneTx_DE/limma_voom_gene_PSD58.R
# Reproduces `volcano_gene_HnrnpM_PSD58_PSDadd.pdf`. Instead of re-running voom, it
#   reads the saved additive gene DE table (limma topTable output: logFC, P.Value,
#   adj.P.Val), which is exactly what the volcano() function plots.
#
# Strategy (same as Fig_2/Fig_3 series): keep ALL original plot parameters untouched,
#   render at the original size, and let ggsave(scale=) squash the figure to
#   7.75 x 6.83 cm (uniform shrink preserves proportions). Pre-inflate ONLY the axis
#   title/tick fonts + legend so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x
#   larger than nominal; resize it to 7.75 cm wide on the slide (that lands the axis
#   at 8 pt). See Fig_2E_exactsize.R for the born-at-size variant.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES      <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
GENE_ADD <- file.path(RES, "02_GeneTx_DE", "limma_gene_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")   # additive ~group+PSD
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 8                     # the source rendered this volcano at width = 8 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Load additive gene DE table (this IS the volcano's data) ----
res <- read.table(GENE_ADD, sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)

# ---- Volcano (original parameters preserved) ----
res$direction <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0, "UP in HnrnpM",
                 ifelse(res$adj.P.Val < 0.05 & res$logFC < 0, "DOWN in HnrnpM", "NS"))
n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
n_up  <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)   # red points
n_dn  <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)   # blue points

pv <- ggplot(res, aes(x = logFC, y = -log10(P.Value), colour = direction)) +
  geom_point(alpha = 0.4, size = 0.6) +
  scale_colour_manual(values = c("UP in HnrnpM"="#d73027",
                                 "DOWN in HnrnpM"="#2166ac","NS"="grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  labs(title = "limma-voom gene: ~group + PSD",
       subtitle = sprintf("%d tested  ·  %d sig at adj.P<0.05 (%d up; %d down)", nrow(res), n_sig, n_up, n_dn),
       x = "log2 Fold Change", y = "-log10(P)", colour = "") +
  guides(colour = guide_legend(override.aes = list(size = 0.6 * 3))) +   # legend dots 3x the plot-dot size
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    legend.position = "bottom",
    axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text  = element_text(size = axis_render_pt)   # -> 8 pt after the ggsave shrink
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_4A.pdf"), pv,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_4A.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_4A.pdf)
# ============================================================
# Volcano plot of gene-level differential expression between control and HNRNPM
# knockdown, under the additive model (~ group + PSD, publication primary). Each
# point is one gene tested with limma-voom on salmon gene counts (filterByExpr +
# TMM normalisation). The horizontal axis is the log2 fold change (HNRNPM knockdown
# relative to control) and the vertical axis is the negative log10 of the unadjusted
# P value. Red points are genes significantly higher in HNRNPM knockdown and blue
# points are genes significantly lower, both at adjusted P below 0.05; grey points
# are not significant. The horizontal dashed line marks an unadjusted P of 0.05 and
# the two vertical dashed lines mark log2 fold changes of minus one and plus one.
# The subtitle reports the number of genes tested and the number significant at
# adjusted P below 0.05.
