# Fig 3A — BSJ vs FSJ logFC scatter, coloured by mechanistic class.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/01_circRNA_DE/FSJ/limma_voom_FSJ_PSD58.R
# Reproduces the `scatter_BSJ_vs_FSJ_logFC.pdf` panel. Instead of re-running the
#   BSJ + FSJ DE pipeline, it reads the saved per-circRNA cross-classification
#   table (BSJ_FSJ_classification.tsv), which already carries logFC_BSJ, logFC_FSJ,
#   sig_BSJ, and the mechanistic class — the exact inputs the scatter uses.
#
# Strategy (same as Fig_2 series): keep ALL original plot parameters untouched,
#   render at the original size, and let ggsave(scale=) squash the figure to
#   7.75 x 6.83 cm (uniform shrink preserves proportions). Pre-inflate ONLY the axis
#   title/tick fonts AND legend fonts so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x
#   larger than 7.75 x 6.83 cm; resize it to 7.75 cm wide on the slide (that also
#   lands the axis/legend at 8 pt). See Fig_2E_exactsize.R for the born-at-size variant.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES      <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
CLASS_F  <- file.path(RES, "01_circRNA_DE", "FSJ", "BSJ_FSJ_classification.tsv")   # primary ~group cross-classification
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 9                     # the source rendered this scatter at width = 9 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Load the cross-classification table (this IS the scatter's data) ----
cross <- read.table(CLASS_F, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
cross$sig_BSJ <- as.logical(cross$sig_BSJ)

# ---- BSJ vs FSJ logFC scatter (original parameters preserved) ----
plot_df <- cross[!is.na(cross$logFC_BSJ) & !is.na(cross$logFC_FSJ), ]
plot_df$facet_cat <- factor(ifelse(plot_df$sig_BSJ, plot_df$class, "BSJ_NS"),
                            levels = c("independent_backsplicing","co_regulated","opposite_regulation",
                                       "FSJ_only_NS_BSJ","neither","BSJ_NS"))
pal2 <- c("independent_backsplicing"="#1f77b4","co_regulated"="#2ca02c",
          "opposite_regulation"="#d62728","FSJ_only_NS_BSJ"="#9467bd",
          "neither"="grey80","BSJ_NS"="grey85")

p_bsj_fsj <- ggplot(plot_df, aes(x = logFC_BSJ, y = logFC_FSJ, colour = facet_cat)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = subset(plot_df, !sig_BSJ), size = 0.4, alpha = 0.25) +
  geom_point(data = subset(plot_df,  sig_BSJ), size = 1.2, alpha = 0.8) +
  scale_colour_manual(values = pal2,
                      labels = c(independent_backsplicing = "Indep",
                                 co_regulated             = "Co-reg",
                                 opposite_regulation      = "Opposite",
                                 FSJ_only_NS_BSJ          = "FSJ_only_NS_BSJ",
                                 neither                  = "neither",
                                 BSJ_NS                   = "BSJ_NS")) +
  labs(title = "BSJ vs FSJ logFC (~group)",
       subtitle = "BSJ-DE circRNAs coloured by mechanistic class",
       x = "BSJ logFC", y = "FSJ logFC", colour = "") +
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    legend.position = "bottom",
    axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text  = element_text(size = axis_render_pt)   # -> 8 pt after the ggsave shrink
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_3A_nonadditive.pdf"), p_bsj_fsj,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_3A_nonadditive.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_3A_nonadditive.pdf)
# ============================================================
# Scatter plot of back-splice junction (BSJ) versus forward-splice junction (FSJ)
# log2 fold change for each circRNA between control and HNRNPM knockdown, under the
# primary group-only model. Each point is one circRNA that could be tested at both
# junctions. The horizontal axis is the BSJ log2 fold change and the vertical axis
# is the FSJ log2 fold change. The dotted lines mark zero on each axis and the
# dashed diagonal marks equal BSJ and FSJ change (y equals x). Grey points are
# circRNAs not significant at the BSJ level; coloured points are BSJ-significant
# circRNAs classified by their mechanism, shown with abbreviated legend labels,
# Indep (independent back-splicing, BSJ changes without a matching FSJ change),
# Co-reg (co-regulated, BSJ and FSJ change together in the same direction), and
# Opposite (opposite regulation, BSJ and FSJ change in opposite directions).
# Points falling near the horizontal (FSJ log2 fold change
# about zero) but away from zero on the horizontal axis represent independent
# back-splicing changes, whereas points along the dashed diagonal represent
# co-regulation of the circular and linear junctions.
