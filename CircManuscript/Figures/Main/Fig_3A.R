# Fig 3A - BSJ vs FSJ logFC scatter, coloured by mechanistic class.
#          ADDITIVE MODEL VERSION (~ group + PSD on BOTH axes).
#
# Additive counterpart of Fig_3A_nonadditive.R (same folder), which uses the
#   ~group model on both axes. This script is identical in structure and in all
#   plot parameters; only the input table, one extra class level, the panel
#   title, the FSJ-untestable handling and the output filenames differ.
#
# TWO OUTPUTS, same data, differing only in how the 4 FSJ-untestable
# ("extreme_linear_collapse") BSJ-DE circRNAs are handled:
#   Fig_3A_v1.pdf - those 4 are SHOWN in a dedicated band below the panel,
#                   with a Collapse legend entry and a subtitle note.
#   Fig_3A_v2.pdf - those 4 are DROPPED: no band, no Collapse legend entry.
#                   The subtitle still records the exclusion so the panel does
#                   not imply 52 is the whole BSJ-DE set.
#
# Input: BSJ_FSJ_classification_additive.tsv, written by
#   .../FinalAnalysis/SCRIPTS/01_circRNA_DE/FSJ/classify_BSJ_FSJ_additive.py
#   from limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv (BSJ axis, ~group + PSD) and
#        limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv (FSJ axis, ~group + PSD).
#   It already carries logFC_BSJ, logFC_FSJ, sig_BSJ and the mechanistic class,
#   so no DE re-run is needed. FSJ-sig uses the publication-adopted SOFT rule
#   (raw P < 0.05 AND |logFC_FSJ| >= log2(1.5)), same rule as the nonadditive panel.
#
# FILENAME NOTE
#   This script no longer writes "Fig_3A.pdf" - it writes Fig_3A_v1.pdf and
#   Fig_3A_v2.pdf. Fig_3A_nonadditive.R writes "Fig_3A_nonadditive.pdf" (its
#   three references were retargeted from "Fig_3A.pdf" when this script was
#   added), so no script in this folder writes "Fig_3A.pdf" any more. Any
#   Fig_3A.pdf still present is a stale artefact of an earlier run.
#
# TUNABLES near the top: TITLE_PT / SUBTITLE_PT set the title and subtitle size
#   (currently HALVED from the previous 9 pt / 7 pt), and DOT_SC scales every
#   scatter mark and the legend key glyph (currently 3x the original sizes).
#
# Strategy (same as Fig_2 series / Fig_3A_nonadditive.R): keep ALL original plot
#   parameters untouched, render at the original size, and let ggsave(scale=)
#   squash the figure to 7.75 x 6.83 cm (uniform shrink preserves proportions).
#   Pre-inflate ONLY the font sizes so they land at their target pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x
#   larger than 7.75 x 6.83 cm; resize it to 7.75 cm wide on the slide (that also
#   lands the axis/legend at 8 pt).

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES      <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
CLASS_F  <- file.path(RES, "01_circRNA_DE", "FSJ", "BSJ_FSJ_classification_additive.tsv")  # additive ~group+PSD cross-classification
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 9                     # the source rendered this scatter at width = 9 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Title / subtitle size ----
# HALVED relative to the previous version of this script (which used 9 pt / 7 pt).
TITLE_PT    <- 9 / 2               # -> 4.5 pt final
SUBTITLE_PT <- 7 / 2               # -> 3.5 pt final

# ---- Scatter mark sizes ----
# DOT_SC scales every scatter mark relative to the original Fig_3A_nonadditive.R
# sizes (0.4 non-sig / 1.2 BSJ-DE / 1.6 band). Set DOT_SC <- 1 to restore them.
DOT_SC      <- 3
SIZE_NS     <- 0.4 * DOT_SC        # grey, BSJ-not-significant
SIZE_SIG    <- 1.2 * DOT_SC        # coloured, BSJ-DE
SIZE_BAND   <- 1.6 * DOT_SC        # FSJ-untestable triangles (v1 only)
SIZE_LEGEND <- SIZE_SIG            # legend key glyph, scaled by the same factor

# ---- Load the cross-classification table (this IS the scatter's data) ----
# na.strings: classify_BSJ_FSJ_additive.py is pandas-written, so the 270 FSJ-filter-fail
# circRNAs carry EMPTY logFC_FSJ / rawP_FSJ / adjP_FSJ fields rather than the string "NA".
cross <- read.table(CLASS_F, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, na.strings = c("NA", "", "nan", "NaN"))
cross$sig_BSJ <- as.logical(cross$sig_BSJ)
cross$sig_BSJ[is.na(cross$sig_BSJ)] <- FALSE

# ---- Split: plottable points vs FSJ-untestable BSJ-DE (extreme linear collapse) ----
# 4 additive BSJ-DE circRNAs (PTK2, BIRC6, SFMBT2, DNAH14) failed the FSJ >=2-in->=2
# abundance filter, so they have NO FSJ logFC and cannot sit on the y axis. v1 shows
# them in a band below the data; v2 drops them. Either way they are never in plot_df.
plot_df <- cross[!is.na(cross$logFC_BSJ) & !is.na(cross$logFC_FSJ), ]
band_df <- cross[cross$sig_BSJ & is.na(cross$logFC_FSJ) & !is.na(cross$logFC_BSJ), ]

n_sig_all  <- sum(cross$sig_BSJ)
n_sig_plot <- sum(plot_df$sig_BSJ)
n_band     <- nrow(band_df)

# ---- Band geometry (derived from the data range; axes stay auto-scaled) ----
y_rng  <- range(plot_df$logFC_FSJ)
y_span <- diff(y_rng)
y_sep  <- y_rng[1] - 0.035 * y_span   # separator rule between panel and band
y_band <- y_rng[1] - 0.090 * y_span   # band point row
x_rng  <- range(plot_df$logFC_BSJ)

CLASS_LEVELS <- c("independent_backsplicing","co_regulated","opposite_regulation",
                  "extreme_linear_collapse","FSJ_only_NS_BSJ","neither","BSJ_NS")

plot_df$facet_cat <- factor(ifelse(plot_df$sig_BSJ, plot_df$class, "BSJ_NS"),
                            levels = CLASS_LEVELS)
band_df$facet_cat <- factor("extreme_linear_collapse", levels = CLASS_LEVELS)
band_df$y_band    <- y_band

# pal2: identical to Fig_3A_nonadditive.R plus extreme_linear_collapse. Note that
# #9467bd is already taken by FSJ_only_NS_BSJ in the nonadditive palette, so the
# collapse class uses #ff7f0e rather than the #9467bd used in Figure1_circRNA_DE.R.
pal2 <- c("independent_backsplicing"="#1f77b4","co_regulated"="#2ca02c",
          "opposite_regulation"="#d62728","extreme_linear_collapse"="#ff7f0e",
          "FSJ_only_NS_BSJ"="#9467bd","neither"="grey80","BSJ_NS"="grey85")

CLASS_LABELS <- c(independent_backsplicing = "Indep",
                  co_regulated             = "Co-reg",
                  opposite_regulation      = "Opposite",
                  extreme_linear_collapse  = "Collapse",
                  FSJ_only_NS_BSJ          = "FSJ_only_NS_BSJ",
                  neither                  = "neither",
                  BSJ_NS                   = "BSJ_NS")

# Shared typography, applied identically to both versions so v1 and v2 are
# typographically indistinguishable.
# PARITY NOTE: Fig_3A_nonadditive.R leaves title/subtitle UNinflated, so they
# shrink to ~4 pt in the final figure. They are inflated here so the
# FSJ-untestable note is actually legible. Delete the plot.title /
# plot.subtitle lines to match Fig_3A_nonadditive.R in appearance.
THEME_FIG <- theme(
  legend.position = "bottom",
  axis.title    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
  axis.text     = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
  legend.title  = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
  legend.text   = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
  plot.title    = element_text(size = TITLE_PT * SC),   # -> 4.5 pt after the shrink
  plot.subtitle = element_text(size = SUBTITLE_PT * SC) # -> 3.5 pt after the shrink
)

PANEL_TITLE <- "BSJ vs FSJ logFC (~group + PSD)"

# ============================================================
# v1 - FSJ-untestable cases SHOWN in a band below the panel
# ============================================================
p_v1 <- ggplot(plot_df, aes(x = logFC_BSJ, y = logFC_FSJ, colour = facet_cat)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = subset(plot_df, !sig_BSJ), size = SIZE_NS, alpha = 0.25) +
  geom_point(data = subset(plot_df,  sig_BSJ), size = SIZE_SIG, alpha = 0.8) +
  # --- FSJ-untestable band ---
  geom_hline(yintercept = y_sep, colour = "grey70", linewidth = 0.3) +
  geom_point(data = band_df, aes(x = logFC_BSJ, y = y_band, colour = facet_cat),
             shape = 17, size = SIZE_BAND, alpha = 0.9, inherit.aes = FALSE) +
  annotate("text", x = x_rng[1], y = y_band, hjust = 0, vjust = -0.9,
           colour = "grey35", size = (7 * SC) / .pt,
           label = sprintf("FSJ untestable (n=%d)", n_band)) +
  scale_colour_manual(values = pal2, labels = CLASS_LABELS) +
  labs(title = PANEL_TITLE,
       subtitle = sprintf("BSJ-DE circRNAs coloured by mechanistic class\n%d/%d plotted - %d FSJ-untestable in band below",
                          n_sig_plot, n_sig_all, n_band),
       x = "BSJ logFC", y = "FSJ logFC", colour = "") +
  guides(colour = guide_legend(override.aes = list(size = SIZE_LEGEND))) +
  theme_bw(base_size = 11) +                          # original base size / proportions
  THEME_FIG

ggsave(file.path(OUT_DIR, "Fig_3A_v1.pdf"), p_v1,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)

# ============================================================
# v2 - the 4 extreme_linear_collapse cases DROPPED
# ============================================================
# plot_df already excludes them (they have no FSJ logFC), so "dropping" means:
# remove the band, the separator rule and the annotation, and remove
# extreme_linear_collapse from the levels/palette/labels so no Collapse entry
# appears in the legend. Nothing about the 5,550 plotted points changes, and the
# y axis now auto-scales to the data instead of reaching down to the band.
CLASS_LEVELS_V2 <- setdiff(CLASS_LEVELS, "extreme_linear_collapse")
plot_df_v2 <- plot_df
plot_df_v2$facet_cat <- factor(as.character(plot_df_v2$facet_cat),
                               levels = CLASS_LEVELS_V2)

p_v2 <- ggplot(plot_df_v2, aes(x = logFC_BSJ, y = logFC_FSJ, colour = facet_cat)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = subset(plot_df_v2, !sig_BSJ), size = SIZE_NS, alpha = 0.25) +
  geom_point(data = subset(plot_df_v2,  sig_BSJ), size = SIZE_SIG, alpha = 0.8) +
  scale_colour_manual(values = pal2[CLASS_LEVELS_V2],
                      labels = CLASS_LABELS[CLASS_LEVELS_V2]) +
  labs(title = PANEL_TITLE,
       subtitle = sprintf("BSJ-DE circRNAs coloured by mechanistic class\n%d/%d plotted - %d FSJ-untestable excluded",
                          n_sig_plot, n_sig_all, n_band),
       x = "BSJ logFC", y = "FSJ logFC", colour = "") +
  guides(colour = guide_legend(override.aes = list(size = SIZE_LEGEND))) +
  theme_bw(base_size = 11) +                          # original base size / proportions
  THEME_FIG

ggsave(file.path(OUT_DIR, "Fig_3A_v2.pdf"), p_v2,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)

cat(sprintf("scale=%.3f; axis/legend rendered %.1f pt -> %d pt final\n",
            SC, axis_render_pt, AXIS_PT))
cat(sprintf("Wrote: %s\n", file.path(OUT_DIR, "Fig_3A_v1.pdf")))
cat(sprintf("Wrote: %s\n", file.path(OUT_DIR, "Fig_3A_v2.pdf")))
cat(sprintf("v1: %d points plotted (%d BSJ-DE) + %d FSJ-untestable BSJ-DE in band = %d BSJ-DE total\n",
            nrow(plot_df), n_sig_plot, n_band, n_sig_all))
cat(sprintf("v2: %d points plotted (%d BSJ-DE); %d FSJ-untestable BSJ-DE dropped\n",
            nrow(plot_df_v2), n_sig_plot, n_band))
cat("Class breakdown among BSJ-DE plotted points (identical in v1 and v2):\n")
print(table(droplevels(plot_df$facet_cat[plot_df$sig_BSJ])))
cat("FSJ-untestable BSJ-DE circRNAs (in v1 band, absent from v2):\n")
print(band_df[, c("circRNA","gene_symbol","logFC_BSJ","adjP_BSJ")], row.names = FALSE)

# ============================================================
# FIGURE LEGEND (for Fig_3A_v1.pdf)
# ============================================================
# Scatter plot of back-splice junction (BSJ) versus forward-splice junction (FSJ)
# log2 fold change for each circRNA between control and HNRNPM knockdown, under the
# additive model that includes the days-post-transduction covariate (~ group + PSD)
# at both junctions. Each point is one circRNA that could be tested at both
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
# co-regulation of the circular and linear junctions. Below the grey rule, the
# band labelled FSJ untestable shows the BSJ-significant circRNAs of class
# Collapse (extreme linear collapse), which have essentially no linear flux at the
# back-splice junction in any sample and therefore fail the FSJ abundance filter;
# they have no FSJ log2 fold change and their horizontal position alone is
# meaningful.
#
# ============================================================
# FIGURE LEGEND (for Fig_3A_v2.pdf)
# ============================================================
# As for Fig_3A_v1.pdf, except that the BSJ-significant circRNAs of class
# extreme linear collapse are omitted. These circRNAs have essentially no linear
# flux at the back-splice junction in any sample, fail the FSJ abundance filter
# and therefore have no FSJ log2 fold change, so they cannot be placed on the
# vertical axis. The panel shows the BSJ-significant circRNAs for which both
# junctions could be tested; the count of omitted circRNAs is stated in the
# panel subtitle.

