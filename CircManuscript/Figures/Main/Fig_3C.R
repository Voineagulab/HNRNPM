# Fig 3C — The 36 independent-back-splicing circRNAs: BSJ vs FSJ CPM box plots.
# Adapted from the `new_36_independent_circ.png` block of
#   /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/Manuscript_R/SCRIPTS/render_new_figures.R,
#   using the FinalAnalysis inputs directly. The source drew ALL 56 BSJ-DE circRNAs
#   with red boxes around the 36 independent ones; THIS version keeps ONLY those 36
#   (class == independent_backsplicing under the additive ~group+PSD soft rule) and
#   drops the other 20 panels, so the grid is a clean 6 x 6.
#
# Strategy (same as Fig_2/Fig_3A-B): keep ALL original plot parameters (base_size 8,
#   box/jitter sizes, strip text) untouched, render at original size, and let
#   ggsave(scale=) squash the figure to 16 x 14.2 cm (uniform shrink preserves
#   proportions). Pre-inflate ONLY the axis title/tick fonts + legend so they land
#   at 8 pt; the facet strip (gene labels) stays at the original proportion.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x
#   larger than 16 x 14.2 cm; resize it to 16 cm wide on the slide (that lands the
#   axis/legend at 8 pt). See Fig_2E_exactsize.R for the born-at-size variant.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
C01      <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE")
CFSJ     <- file.path(C01, "FSJ")
CLASS_ADD <- file.path(CFSJ, "BSJ_FSJ_classification_additive.tsv")   # additive soft-rule classification
CPM_BSJ  <- file.path(C01,  "cpm_BSJ_HnrnpM_PSD58.tsv")
CPM_FSJ  <- file.path(CFSJ, "cpm_FSJ_HnrnpM_PSD58.tsv")
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 16; H_CM <- 14.2        # nominal final PDF size on the slide
ORIG_W_IN <- 12                    # the source rendered the panel grid at width = 12 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Select the 36 independent-back-splicing circRNAs (additive) ----
rd <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "",
                             check.names = FALSE, stringsAsFactors = FALSE)
cls <- rd(CLASS_ADD)
indep <- cls[cls$class == "independent_backsplicing", ]
indep <- indep[order(indep$adjP_BSJ), ]              # most significant first
cat("independent-back-splicing circRNAs:", nrow(indep), "\n")

# Unique gene-name panel labels (disambiguate duplicate host genes with #index)
gl  <- indep$gene_symbol
dup <- ave(seq_along(gl), gl, FUN = length) > 1
idx <- ave(seq_along(gl), gl, FUN = seq_along)
indep$panel_title <- ifelse(dup, paste0(gl, "#", idx), gl)

# ---- Long-format BSJ + FSJ CPM for those 36 ----
bsj_cpm <- rd(CPM_BSJ); fsj_cpm <- rd(CPM_FSJ)
neg_cols <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
hm_cols  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
circs    <- indep$circRNA
mk_long <- function(cpm_df, junction) {
  d <- cpm_df[cpm_df$circRNA %in% circs, ]
  do.call(rbind, lapply(c(neg_cols, hm_cols), function(s)
    data.frame(circRNA = d$circRNA, cpm = d[[s]],
               group = ifelse(s %in% hm_cols, "gHNRNPM", "gCTRL"),
               junction = junction, stringsAsFactors = FALSE)))
}
pdat <- rbind(mk_long(bsj_cpm, "BSJ"), mk_long(fsj_cpm, "FSJ"))
pdat$log2cpm <- log2(pdat$cpm + 1)
pdat$jg <- factor(paste(pdat$junction, pdat$group, sep = "-"),
                  levels = c("BSJ-gCTRL","BSJ-gHNRNPM","FSJ-gCTRL","FSJ-gHNRNPM"))
pdat$panel_title <- factor(setNames(indep$panel_title, indep$circRNA)[pdat$circRNA],
                           levels = indep$panel_title)

pal <- c("BSJ-gCTRL"="#4C78A8","BSJ-gHNRNPM"="#E45756","FSJ-gCTRL"="#9ec3e3","FSJ-gHNRNPM"="#f5a8a8")

p36 <- ggplot(pdat, aes(jg, log2cpm, fill = jg)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, colour = "grey20") +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.8, colour = "grey20") +
  scale_fill_manual(values = pal, name = NULL) +
  facet_wrap(~ panel_title, ncol = 6, scales = "free_y") +
  labs(title = sprintf("Independent back-splicing circRNAs (additive ~ group + PSD): BSJ vs FSJ CPM  (n = %d)", nrow(indep)),
       x = NULL, y = "log2(CPM + 1)") +
  theme_bw(base_size = 8) +                                     # original base size / proportions
  theme(
    strip.text      = element_text(size = 6.6 * 1.5),           # gene labels: 1.5x original
    axis.text.x     = element_text(angle = 45, hjust = 1, size = axis_render_pt),  # -> 8 pt after shrink
    axis.text.y     = element_text(size = axis_render_pt),                          # -> 8 pt after shrink
    axis.title      = element_text(size = axis_render_pt),                          # -> 8 pt after shrink
    legend.text     = element_text(size = axis_render_pt),                          # -> 8 pt after shrink
    legend.position = "bottom",
    panel.spacing   = grid::unit(0.4, "lines")
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_3C.pdf"), p36,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_3C.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_3C.pdf)
# ============================================================
# Per-circRNA back-splice junction (BSJ) and forward-splice junction (FSJ)
# expression for the 36 independent-back-splicing circRNAs, the BSJ-DE circRNAs
# whose back-splice product changes under HNRNPM knockdown without a matching change
# in the linear (FSJ) junction (additive ~ group + PSD model; FSJ classed as
# unchanged under the soft rule of raw P >= 0.05 or |FSJ log2 fold change| <
# log2(1.5)). Each small panel is one circRNA, labelled by its host gene (duplicate
# host genes disambiguated with a #index). Within each panel the four boxes show
# log2(CPM + 1) for BSJ in control (gCTRL) and HNRNPM knockdown (gHNRNPM) and FSJ in
# control and knockdown, with individual samples overlaid as points. Boxes show the
# median and interquartile range. In this class the BSJ boxes shift between control
# and knockdown while the FSJ boxes at the same junction stay comparable, the
# signature of a change in back-splicing that is independent of the host gene's
# linear splicing.
