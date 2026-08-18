# Fig 3B — BSJ x FSJ mechanistic class barplot (additive ~group+PSD model only).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/05_Summary_Plots/BSJ_FSJ_classification_barplot_PSD58.R
# Retains only the LEFT facet of that figure — the additive-model (~ group + PSD,
#   publication primary) classification of BSJ-DE circRNAs into 4 mechanistic
#   classes. The primary (~group) comparator facet is dropped.
#
# Strategy (same as Fig_2/Fig_3A): keep ALL original plot parameters (base_size 11,
#   bar width, bar-label size) untouched so the plot looks IDENTICAL to the source
#   panel. Render at the original per-panel size and let ggsave(scale=) squash the
#   figure to 7.75 x 6.83 cm (uniform shrink preserves proportions). Pre-inflate
#   ONLY the axis title/tick fonts so they land at 8 pt (fill legend is hidden).
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x
#   larger than nominal; resize it to 7.75 cm wide on the slide (that lands the axis
#   at 8 pt). See Fig_2E_exactsize.R for the born-at-size variant.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
FSJ_DIR <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ")
ADD_TBL <- file.path(FSJ_DIR, "BSJ_FSJ_classification_additive.tsv")   # additive ~group+PSD
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 5.5                   # each facet of the source 2-panel figure (11 in wide) is ~5.5 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Load additive classification, restrict to BSJ-DE circRNAs ----
to_bool <- function(x) x == "True" | x == TRUE
add <- read.table(ADD_TBL, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
add$sig_BSJ <- to_bool(add$sig_BSJ)
add_sig <- add[add$sig_BSJ, ]

# ---- Tabulate the 4 mechanistic classes ----
class_levels <- c("independent_backsplicing","co_regulated","opposite_regulation","extreme_linear_collapse")
n_by_class <- sapply(class_levels, function(k) sum(add_sig$class == k))
total <- sum(n_by_class)
df <- data.frame(
  class = factor(class_levels, levels = class_levels),
  n     = as.integer(n_by_class),
  pct   = round(100 * n_by_class / total, 1)
)

pal <- c("independent_backsplicing"="#1f77b4","co_regulated"="#2ca02c",
         "opposite_regulation"="#d62728","extreme_linear_collapse"="#9467bd")

p <- ggplot(df, aes(x = class, y = pct, fill = class)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)), vjust = -0.4, size = 3) +  # original bar-label size
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = c(independent_backsplicing = "Indep",
                              co_regulated             = "Co-reg",
                              opposite_regulation      = "Opposite",
                              extreme_linear_collapse  = "Collapse")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = sprintf("BSJ x FSJ mechanistic class (additive ~group+PSD, n=%d)", total),
       subtitle = "Independent back-splicing dominates the HNRNPM-BSJ-DE circRNAs.",
       x = NULL, y = "% of BSJ-DE circRNAs") +
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5, size = axis_render_pt),  # horizontal, -> 8 pt after shrink
    axis.text.y  = element_text(size = axis_render_pt),                          # -> 8 pt after shrink
    axis.title   = element_text(size = axis_render_pt),                          # -> 8 pt after shrink
    plot.subtitle = element_text(size = 5 * SC)                                  # -> ~5 pt after shrink (smaller)
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_3B.pdf"), p,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_3B.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_3B.pdf)
# ============================================================
# Mechanistic classification of the HNRNPM BSJ-differentially-expressed circRNAs
# under the additive model (~ group + PSD, publication primary). Each BSJ-DE circRNA
# is assigned to one of four classes by comparing its back-splice junction (BSJ) and
# forward-splice junction (FSJ) differential expression, independent back-splicing
# (BSJ changes without a significant FSJ change), co-regulated (BSJ and FSJ change
# together in the same direction), opposite regulation (BSJ and FSJ change in
# opposite directions), and extreme linear collapse (the linear FSJ is essentially
# absent so FSJ differential expression cannot be tested). Bars show the percentage
# of BSJ-DE circRNAs in each class, with the count and percentage printed above each
# bar. The horizontal axis is the mechanistic class and the vertical axis is the
# percentage of BSJ-DE circRNAs. The dominance of the independent back-splicing
# class indicates that most HNRNPM-responsive circRNA changes occur at the
# back-splice junction independently of the host gene's linear splicing.
