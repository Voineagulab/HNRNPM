# Fig 2D — CLR distribution by bin (barplot), bottom-right panel of the CLR figure.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/03_CLR/global_shift_plots_PSD58.R
# Retains only the 4-bin CLR histogram panel — the % of detected circRNAs per CLR
#   bin (<0.1 / 0.1-0.5 / 0.5-0.9 / >0.9), dodged by group. The scatter, density,
#   and by-sample density panels are dropped.
#
# Strategy (same as Fig_2A-2D): keep ALL original plot parameters (base_size 11,
#   bar widths, bar-label size) untouched so the plot looks IDENTICAL to the source.
#   Render at the original size and let ggsave(scale=) squash the figure down to
#   7.75 x 6.83 cm — a uniform shrink preserves every proportion. Pre-inflate ONLY the
#   axis title/tick fonts AND legend fonts so they land at 8 pt in the output.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(tidyr); library(dplyr); library(tibble)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES     <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
CLR_F   <- file.path(RES, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83        # final PDF size on the slide
ORIG_W_IN <- 9                     # the source rendered this panel at width = 9 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Load CLR matrix ----
dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE)
meta_cols <- c("circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id")
mat <- as.matrix(dat[, setdiff(colnames(dat), meta_cols)])
rownames(mat) <- dat$circRNA

# ---- Long format, detected only (CLR > 0), group relabelled gCTRL/gHNRNPM ----
long <- as.data.frame(mat) %>% tibble::rownames_to_column("circRNA") %>%
  pivot_longer(-circRNA, names_to = "sample", values_to = "CLR") %>%
  mutate(group = factor(ifelse(grepl("gHnrnpM", sample), "gHNRNPM", "gCTRL"),
                        levels = c("gCTRL","gHNRNPM"))) %>%
  filter(CLR > 0)

group_pal <- c(gCTRL = "#4C78A8", gHNRNPM = "#E45756")

# ---- 4-bin distribution (% of detected per group) ----
bin_breaks <- c(0, 0.1, 0.5, 0.9, 1.001)
bin_labels <- c("<0.1", "0.1-0.5", "0.5-0.9", ">0.9")
bin_df <- long %>%
  mutate(bin = cut(CLR, breaks = bin_breaks, labels = bin_labels,
                   include.lowest = TRUE, right = FALSE)) %>%
  dplyr::count(group, bin, name = "n") %>%           # namespace-qualified to avoid masking
  dplyr::group_by(group) %>% dplyr::mutate(pct = 100 * n / sum(n)) %>% dplyr::ungroup()

pe <- ggplot(bin_df, aes(x = bin, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "white") +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 2.6) +   # original bar-label size
  scale_fill_manual(values = group_pal) +
  labs(title = "CLR distribution by bin (% of detected, per group)",
       x = "CLR bin", y = "% detected per group", fill = "Group") +
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    legend.position = "right",
    axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text  = element_text(size = axis_render_pt)   # -> 8 pt after the ggsave shrink
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_2D.pdf"), pe,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_2D.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_2D.pdf)
# ============================================================
# Distribution of circRNAs across four circular-to-linear ratio (CLR) bins,
# compared between control and HNRNPM knockdown. The CLR of a circRNA in a sample
# is 2 times its back-splice junction count divided by the sum of 2 times the
# back-splice junction count and the forward-splice junction count, ranging from 0
# (fully linear) to 1 (fully circular). Each detected circRNA-sample observation
# (CLR above 0) is assigned to one of four bins (below 0.1, 0.1 to 0.5, 0.5 to 0.9,
# and above 0.9). Bars show the percentage of detected observations in each bin,
# computed separately within each group and dodged side by side, control (gCTRL) in
# blue and HNRNPM knockdown (gHNRNPM) in red, with the exact percentage printed
# above each bar. The horizontal axis is the CLR bin and the vertical axis is the
# percentage of detected observations per group. A shift of mass out of the lowest
# bin into the higher bins under HNRNPM knockdown reflects the transcriptome-wide
# increase in circular fraction.
