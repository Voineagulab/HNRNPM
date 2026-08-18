# Fig 2B — CLR density by sample (per-sample curves, gCTRL dashed / gHNRNPM solid).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/03_CLR/global_shift_plots_PSD58.R
# Retains only the by-SAMPLE density panel (middle-right of the combined figure) —
#   the distribution of circular-to-linear ratio (CLR) for detected circRNAs, one
#   curve per sample, coloured by sample and dashed/solid by group, on a log10 x
#   axis. The scatter, by-group density, and histogram panels are dropped.
#
# Strategy (same as Fig_2A/2B): keep ALL original plot parameters (base_size 11,
#   density line width, palettes) untouched so the plot looks IDENTICAL to the
#   source. Render at the original size and let ggsave(scale=) squash the whole
#   figure down to 7.75 x 6.83 cm — a uniform shrink preserves every proportion.
#   Pre-inflate ONLY the axis title/tick fonts AND the legend fonts by the same
#   factor so they land at exactly 8 pt in the output.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(tidyr); library(dplyr); library(tibble); library(RColorBrewer)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES     <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
CLR_F   <- file.path(RES, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # final PDF size on the slide
ORIG_W_IN <- 9                     # the source rendered this density panel at width = 9 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Sample groups (order: controls first, then knockdown) ----
HM  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
NEG <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
sample_order <- c(NEG, HM)

# ---- Load CLR matrix ----
dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE)
meta_cols <- c("circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id")
mat <- as.matrix(dat[, setdiff(colnames(dat), meta_cols)])
rownames(mat) <- dat$circRNA

# ---- Long format, detected only (CLR > 0) ----
long <- as.data.frame(mat) %>% tibble::rownames_to_column("circRNA") %>%
  pivot_longer(-circRNA, names_to = "sample", values_to = "CLR") %>%
  mutate(group  = factor(ifelse(grepl("gHnrnpM", sample), "gHNRNPM", "gCTRL"),
                         levels = c("gCTRL","gHNRNPM")),
         sample = factor(sample, levels = sample_order)) %>%
  filter(CLR > 0)

# Per-sample colour palette (Blues for controls, Reds for knockdown) — keyed by the
# real sample names so the mapping matches the data.
sample_pal <- c(setNames(brewer.pal(9, "Blues")[5:9][seq_along(NEG)], NEG),
                setNames(brewer.pal(7, "Reds")[4:7][seq_along(HM)],   HM))

# Relabel the sample legend entries gNEG4 -> gCTRL, gHnrnpM -> gHNRNPM (display only)
relabel <- function(x) unname(c("9_gNEG4_PSD5_S24" = "gCTRL_PSD5_rep1", "12_gNEG4_PSD5_S4" = "gCTRL_PSD5_rep2", "11_gNEG4_PSD8_S3" = "gCTRL_PSD8_rep1", "13_gNEG4_PSD8_S5" = "gCTRL_PSD8_rep2", "5_gHnrnpM_PSD5_S20" = "gHNRNPM_PSD5_rep1", "6_gHnrnpM_PSD5_S21" = "gHNRNPM_PSD5_rep2", "7_gHnrnpM_PSD8_S22" = "gHNRNPM_PSD8")[as.character(x)])

pd <- ggplot(long, aes(x = CLR, colour = sample, linetype = group, group = sample)) +
  geom_density(linewidth = 0.7) +                                # original line width
  scale_colour_manual(values = sample_pal, labels = relabel, name = "Sample") +
  scale_linetype_manual(values = c(gCTRL = "dashed", gHNRNPM = "solid")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), limits = c(0.001, 1)) +
  labs(title = "CLR density by sample (log-x)",
       subtitle = "Blues = gCTRL (dashed); Reds = gHNRNPM (solid)",
       x = "CLR (log scale)", y = "Density", colour = "Sample", linetype = "Group") +
  guides(
    colour   = guide_legend(position = "bottom", nrow = 4, byrow = FALSE,   # 4 per column: controls left, KD right
                            theme = theme(legend.text = element_text(size = 5 * SC))),          # Sample entries -> ~5 pt final
    linetype = guide_legend(position = "right",
                            theme = theme(legend.text = element_text(size = axis_render_pt)))   # Group entries  -> 8 pt final
  ) +
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    legend.position = "right",
    axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text  = element_text(size = axis_render_pt)   # -> 8 pt after the ggsave shrink
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_2B.pdf"), pd,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_2B.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_2B.pdf)
# ============================================================
# Per-sample density distribution of the circular-to-linear ratio (CLR) across
# circRNAs. The CLR of a circRNA in a sample is 2 times its back-splice junction
# count divided by the sum of 2 times the back-splice junction count and the
# forward-splice junction count, ranging from 0 (fully linear) to 1 (fully
# circular). Each observation is one circRNA in one sample, restricted to detected
# entries with a CLR above 0. The horizontal axis shows CLR on a base 10
# logarithmic scale and the vertical axis shows kernel density. One curve is drawn
# per sample, coloured in blues for the four control samples (gCTRL) and reds for
# the three HNRNPM knockdown samples (gHNRNPM); line type additionally encodes
# group, dashed for control and solid for knockdown. The consistent rightward
# position of the knockdown (red, solid) curves relative to the control (blue,
# dashed) curves shows that the shift toward higher circular fraction under HNRNPM
# knockdown is reproducible across individual replicates, not driven by a single
# sample.
