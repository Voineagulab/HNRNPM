# Fig 4E — Transcript-level volcano (additive ~ group + PSD).
# Adapted from the "transcript volcano rebuilt from data" block of
#   /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/Manuscript_R/SCRIPTS/render_new_figures.R
#   (which rebuilds panel_transcript_volcano.png), reading the FinalAnalysis 02
#   additive transcript DE table directly.
#
# Strategy (same as the Fig_2/Fig_3/Fig_4 series): keep ALL original plot parameters
#   (point size/alpha, colours, thresholds, title/subtitle sizes), render at the
#   original size, and let ggsave(scale=) squash the figure to 7.75 x 6.83 cm (uniform
#   shrink preserves proportions). Pre-inflate ONLY the axis title/tick fonts + legend
#   so they land at 8 pt in the output.
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
TX_ADD   <- file.path(RES, "02_GeneTx_DE", "limma_transcript_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")   # additive ~group+PSD
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 8                     # the source rendered this volcano at width = 8 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Load additive transcript DE table (this IS the volcano's data) ----
tx <- read.table(TX_ADD, header = TRUE, sep = "\t", quote = "", check.names = FALSE,
                 stringsAsFactors = FALSE)

# ---- Volcano (original parameters preserved) ----
tx$direction <- ifelse(tx$adj.P.Val < 0.05 & tx$logFC > 0, "UP in HnrnpM",
                ifelse(tx$adj.P.Val < 0.05 & tx$logFC < 0, "DOWN in HnrnpM", "NS"))
n_sig <- sum(tx$adj.P.Val < 0.05, na.rm = TRUE)
n_up  <- sum(tx$direction == "UP in HnrnpM")
n_dn  <- sum(tx$direction == "DOWN in HnrnpM")

ptx <- ggplot(tx, aes(x = logFC, y = -log10(P.Value), colour = direction)) +
  geom_point(alpha = 0.3, size = 0.4) +
  scale_colour_manual(values = c("UP in HnrnpM" = "#d73027",
                                 "DOWN in HnrnpM" = "#2166ac", "NS" = "grey80")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  labs(title = "limma-voom transcript: ~group + PSD",
       subtitle = sprintf("%d tested  ·  %d sig at adj.P<0.05  (UP %d / DOWN %d)",
                          nrow(tx), n_sig, n_up, n_dn),
       x = "log2 Fold Change", y = "-log10(P)", colour = "") +
  guides(colour = guide_legend(override.aes = list(size = 0.4 * 3))) +   # legend dots 3x the plot dots
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    legend.position = "bottom",
    axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text  = element_text(size = axis_render_pt)   # -> 8 pt after the ggsave shrink
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_4E.pdf"), ptx,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_4E.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_4E.pdf)
# ============================================================
# Volcano plot of transcript-level differential expression between control and
# HNRNPM knockdown, under the additive model (~ group + PSD). Each point is one
# transcript tested with limma-voom on salmon transcript counts. The horizontal axis
# is the log2 fold change (HNRNPM knockdown relative to control) and the vertical
# axis is the negative log10 of the unadjusted P value. Red points are transcripts
# significantly higher in HNRNPM knockdown and blue points significantly lower, both
# at adjusted P below 0.05; grey points are not significant. The horizontal dashed
# line marks an unadjusted P of 0.05 and the two vertical dashed lines mark log2 fold
# changes of minus one and plus one. The subtitle reports the number of transcripts
# tested, the number significant at adjusted P below 0.05, and the up/down split of
# the coloured (significant) transcripts.
