# Fig 2Z — CLR per-sample violins (detected CLR, log10) with the per-sample mean CLR.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d2_CLR/d2_CLR_violins.R
# Retains only PLOT 1 of that script — one violin per sample (4 gCTRL, 3 gHNRNPM) of the
#   circular-to-linear ratio (CLR) of detected circRNAs on a log10 axis, with a black
#   segment marking that sample's mean CLR. PLOT 2 (the per-group paired violins) is
#   dropped. The Welch t / Mann-Whitney U P values annotated in the subtitle are read
#   from global_CLR_shift_tests.tsv so the panel stays faithful to the tested numbers.
#
# Strategy (same as Fig_2A-2D): keep ALL original plot parameters (base_size 9, violin
#   scale/trim, line widths, mean-line geometry) untouched so the panel looks IDENTICAL
#   to the source. Render at the original size and let ggsave(scale=) squash the figure
#   down to 7.75 x 6.83 cm — a uniform shrink preserves every proportion. Pre-inflate the
#   fonts by the same factor so they land at exactly 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger
#   than nominal; dragging it onto the slide and setting the size to 7.75 x 6.83 cm
#   applies the 1/SC shrink that lands the text at 8 pt.
# NOTE: the axis, legend and in-plot mean labels are pinned to 8 pt final; the title and
#   subtitle are pinned to AXIS_PT/5 = 1.6 pt, i.e. one fifth of that, so they stay
#   subordinate to the data (the sibling panels sit in the same range — Fig_1C uses 3 and
#   4 pt, Fig_3A 4.5 and 3.5 pt).
# NOTE: at 8 pt a 7.75 cm panel fits only ~6 characters per x-axis slot, so the sample
#   identity is split between the x axis (timepoint + replicate, two lines) and the fill
#   legend (group), as in Fig_2D; anything wider is silently clipped by cairo_pdf. The
#   subtitle is kept wrapped to three short lines for the same reason, which now leaves it
#   comfortably inside the panel at 1.6 pt.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(tidyr); library(dplyr); library(tibble)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES     <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
CLR_F   <- file.path(RES, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
TESTS   <- file.path(RES, "03_CLR", "global_CLR_shift_tests.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 18 / 2.54             # the source rendered this violin panel at width = 18 cm
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend + mean-label font size
TITLE_PT <- AXIS_PT / 5            # desired FINAL title + subtitle size = 1/5 of the above
axis_render_pt  <- AXIS_PT * SC           # pre-inflate so it becomes 8 pt after the 1/SC shrink
title_render_pt <- TITLE_PT * SC          # pre-inflate so it becomes 1.6 pt after the shrink
label_render_mm <- axis_render_pt / .pt   # same 8 pt final, in mm for geom_text

# ---- Sample groups (order: controls first, then knockdown — same order as Fig_2B) ----
HM  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
NEG <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
sample_order <- c(NEG, HM)

# ---- Load CLR matrix ----
dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE)
meta_cols <- c("circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id")
mat <- as.matrix(dat[, setdiff(colnames(dat), meta_cols)])
rownames(mat) <- dat$circRNA

# ---- Read the tests table so the annotated P values stay faithful to the analysis ----
tests <- read.table(TESTS, sep = "\t", header = TRUE, check.names = FALSE)
getp <- function(pat) {
  hit <- tests$p[grepl(pat, tests$test)]
  if (length(hit) == 0) NA_real_ else as.numeric(hit[1])
}
fp <- function(p) {                                  # pretty p-value
  if (is.na(p)) return("NA")
  if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "g", digits = 2)
}
p_welch <- getp("Welch_t")
p_mwu   <- getp("Mann-Whitney")

# ---- Long format; group relabelled gCTRL/gHNRNPM ----
long <- as.data.frame(mat) %>% tibble::rownames_to_column("circRNA") %>%
  pivot_longer(-circRNA, names_to = "sample", values_to = "CLR") %>%
  mutate(group  = factor(ifelse(grepl("gHnrnpM", sample), "gHNRNPM", "gCTRL"),
                         levels = c("gCTRL","gHNRNPM")),
         sample = factor(sample, levels = sample_order))

# Violins: detected entries only (CLR = 0 cannot be shown on a log axis)
long_det <- dplyr::filter(long, CLR > 0)

# Mean line: arithmetic mean over ALL 5,820 circRNAs (including the zeros) — that is the
# exact quantity compared by the per-sample Welch t / Mann-Whitney U tests annotated below.
samp_mean <- long %>%
  dplyr::group_by(sample, group) %>%
  dplyr::summarise(mean_CLR = mean(CLR), .groups = "drop") %>%
  dplyr::mutate(x = as.integer(sample))

group_pal  <- c(gCTRL = "#4C78A8", gHNRNPM = "#E45756")
log_breaks <- c(0.001, 0.01, 0.1, 1)

# Relabel the sample axis entries (display only). Same tokens as the Fig_2B sample legend,
# minus the group prefix (the fill legend carries the group) and wrapped onto two lines so
# seven labels fit unrotated at 8 pt.
relabel <- function(x) unname(c("9_gNEG4_PSD5_S24" = "PSD5\nrep1", "12_gNEG4_PSD5_S4" = "PSD5\nrep2", "11_gNEG4_PSD8_S3" = "PSD8\nrep1", "13_gNEG4_PSD8_S5" = "PSD8\nrep2", "5_gHnrnpM_PSD5_S20" = "PSD5\nrep1", "6_gHnrnpM_PSD5_S21" = "PSD5\nrep2", "7_gHnrnpM_PSD8_S22" = "PSD8")[as.character(x)])

pz <- ggplot(long_det, aes(x = sample, y = CLR, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey30",
              linewidth = 0.3, alpha = 0.6) +                      # original violin params
  geom_segment(data = samp_mean,
               aes(x = x - 0.45, xend = x + 0.45, y = mean_CLR, yend = mean_CLR),
               inherit.aes = FALSE, linewidth = 0.8, colour = "black") +
  geom_text(data = samp_mean,
            aes(x = x, y = mean_CLR, label = sprintf("%.3f", mean_CLR)),
            inherit.aes = FALSE, vjust = -0.6, size = label_render_mm) +   # -> 8 pt final
  scale_fill_manual(values = group_pal, name = "Group") +
  scale_x_discrete(labels = relabel) +
  scale_y_log10(breaks = log_breaks, labels = log_breaks) +
  labs(title = "CLR per sample (log10 scale)",
       subtitle = sprintf("Violins: detected CLR (>0); bar = mean\nWelch t P = %s (3 KD vs 4 gCTRL)\nMann-Whitney U P = %s",
                          fp(p_welch), fp(p_mwu)),
       x = NULL, y = "CLR (log10 scale)") +
  theme_bw(base_size = 9) +                           # original base size / proportions
  theme(
    legend.position = "bottom",
    plot.title    = element_text(size = title_render_pt), # -> 1.6 pt after the shrink (8/5)
    plot.subtitle = element_text(size = title_render_pt), # -> 1.6 pt after the shrink (8/5)
    axis.title    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text     = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title  = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text   = element_text(size = axis_render_pt)   # -> 8 pt after the ggsave shrink
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_2Z.pdf"), pz,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis rendered %.1f pt -> %d pt final; title/subtitle %.1f pt -> %.1f pt final)\n",
            file.path(OUT_DIR, "Fig_2Z.pdf"), SC, axis_render_pt, AXIS_PT,
            title_render_pt, TITLE_PT))

# ============================================================
# FIGURE LEGEND (for Fig_2Z.pdf)
# ============================================================
# Per-sample distribution of the circular-to-linear ratio (CLR) across circRNAs. The CLR
# of a circRNA in a sample is 2 times its back-splice junction count divided by the sum of
# 2 times the back-splice junction count and the forward-splice junction count, ranging
# from 0 (fully linear) to 1 (fully circular). One violin is drawn per sample, four
# control samples (gCTRL, blue) followed by three HNRNPM knockdown samples (gHNRNPM, red),
# with violins scaled to equal width; the horizontal axis gives each sample's
# post-selection day and replicate number and the fill colour gives its group. The
# vertical axis shows CLR on a base 10 logarithmic
# scale; because a CLR of zero cannot be placed on a logarithmic axis, the violins show
# only detected entries (CLR above 0). The black horizontal bar within each violin marks
# that sample's arithmetic mean CLR computed across all 5,820 circRNAs, including
# undetected entries scored as zero, and the printed number gives its value; this is the
# quantity compared between groups by the tests reported in the subtitle. The subtitle
# reports the one-sided Welch t-test and Mann-Whitney U P values for higher mean CLR in
# knockdown than in control, taking the sample as the unit of replication (three knockdown
# versus four control). The consistently higher position of the knockdown means, with no
# overlap between the two groups' ranges, shows a reproducible transcriptome-wide shift
# toward higher circular fraction upon HNRNPM knockdown.
