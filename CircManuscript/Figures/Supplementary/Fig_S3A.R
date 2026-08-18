# Fig S3A — FSJ differential-expression volcano, HnrnpM vs NEG (additive ~group + PSD).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/01_circRNA_DE/FSJ/limma_voom_FSJ_PSD58.R
# Reproduces `volcano_FSJ_HnrnpM_PSD58_PSDadd.pdf` (the additive-model FSJ volcano, the
#   `volcano(res_add, "limma-voom FSJ: ~group + PSD")` call). Rather than re-running the
#   limma-voom pipeline, it reads the saved additive FSJ DE table directly.
#
# Strategy (same as the Fig_2/Fig_3/Fig_4/Fig_5 series): keep ALL original plot parameters
#   (point size/alpha, dashed threshold lines, title/subtitle, bottom legend), render at the
#   original size, and let ggsave(scale=) squash the figure to 7.75 x 6.83 cm (uniform shrink
#   preserves proportions). Pre-inflate ONLY the axis title/tick fonts + legend so they land
#   at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger than
#   nominal; resize it to 7.75 cm wide on the slide. See Fig_2E_exactsize.R.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RESFSJ  <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ")
DE_TSV  <- file.path(RESFSJ, "limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")   # additive (~group + PSD)
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Supplementary")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 8                     # the source rendered the volcano at width = 8 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- load the additive FSJ DE table ----
res <- read.table(DE_TSV, sep = "\t", header = TRUE, check.names = FALSE, quote = "")

# ---- direction + counts (identical to source volcano()) ----
res$direction <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0, "UP in HnrnpM",
                 ifelse(res$adj.P.Val < 0.05 & res$logFC < 0, "DOWN in HnrnpM", "NS"))
n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
n_up  <- sum(res$direction == "UP in HnrnpM")
n_dn  <- sum(res$direction == "DOWN in HnrnpM")

# ---- build the volcano (original parameters preserved) ----
pg <- ggplot(res, aes(x = logFC, y = -log10(P.Value), colour = direction)) +
  geom_point(alpha = 0.5, size = 0.9) +
  scale_colour_manual(values = c("UP in HnrnpM"="#d73027",
                                  "DOWN in HnrnpM"="#2166ac","NS"="grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  labs(title = "limma-voom FSJ: ~group + PSD",
       subtitle = sprintf("%d tested  ·  %d sig (UP %d / DOWN %d) at adj.P<0.05",
                          nrow(res), n_sig, n_up, n_dn),
       x = "log2 Fold Change", y = "-log10(P)", colour = "") +
  guides(colour = guide_legend(override.aes = list(size = 0.9 * 3))) +   # legend dots 3x the plot dots
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.title  = element_text(size = axis_render_pt),   # -> 8 pt after shrink
        axis.text   = element_text(size = axis_render_pt),   # -> 8 pt after shrink
        legend.text = element_text(size = axis_render_pt))   # -> 8 pt after shrink

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_S3A.pdf"), pg,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_S3A.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_S3A.pdf)
# ============================================================
# Volcano plot of forward-splice-junction (FSJ) differential expression between HNRNPM
# knockdown and control, from the additive limma-voom model (~ group + PSD). Each point
# is one circRNA locus' FSJ; the horizontal axis is the log2 fold change (knockdown vs
# control) and the vertical axis is -log10 of the raw P value. Points are coloured by
# significance and direction at adjusted P < 0.05 (red up in knockdown, blue down in
# knockdown, grey not significant). The horizontal dashed line marks the nominal
# P = 0.05 level and the vertical dashed lines mark log2 fold changes of -1 and +1. The
# subtitle reports the number of FSJs tested and the number significant (up and down) at
# adjusted P < 0.05.
