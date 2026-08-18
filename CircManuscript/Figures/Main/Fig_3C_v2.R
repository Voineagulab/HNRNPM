# Fig 3C v2 — IDENTICAL to Fig_3C_v1.R except that every individual panel now carries a
#             FULLY INDEPENDENT y-axis (each circRNA's BSJ and FSJ panels are scaled on
#             their own log2(CPM+1) range), instead of v1's per-junction-row shared y.
#
# Exports ONE file: Figures/Main/Fig_3C_v2.pdf
#
# WHY THE PLOT ASSEMBLY DIFFERS FROM v1
#   v1 draws each 12-circRNA block as one facet_grid(junction ~ gene). facet_grid's
#   scales = "free_y" frees y only ACROSS ROWS, so all 12 columns in a junction row are
#   forced onto a shared y-scale. A single ggplot facet cannot give every panel its own
#   y-axis without ggh4x (not installed here). So v2 renders each circRNA as its OWN
#   single-column facet_grid(junction ~ gene, scales = "free_y") — with one column,
#   free_y frees the two junction rows against each other, i.e. an independent y-axis for
#   that circRNA's BSJ and FSJ panel — and lays 12 of them side by side per block with
#   patchwork. Across circRNAs the plots are separate, so all panels are independent.
#
# The 36-circRNA additive independent_backsplicing set (ordered by adjP_BSJ, #index
#   de-duplication), 2-colour group fill, PSD point shape, grey strips, 8 pt fonts,
#   single bottom legend, 16 cm width and cairo_pdf are unchanged from Fig_3C_v1.R.
#   LAYOUT here: four stacked blocks of nine columns (36 = 4 x 9); y-axis tick labels
#   sit on the LEFT of each panel (title once at the far left); the BSJ/FSJ grey strips
#   sit on the RIGHT of the last panel in each row; the x-axis group tick labels are
#   removed as redundant with the fill legend. Rendered at Fig_3C.pdf's page size (see
#   the sizing block) so the fonts match Fig_3C.pdf once resized to 16 x 14.2 cm.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})


# --- Path definition ---
ROOT      <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
C01       <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE")
CFSJ      <- file.path(C01, "FSJ")
CLASS_ADD <- file.path(CFSJ, "BSJ_FSJ_classification_additive.tsv")   # additive soft-rule classification
CPM_BSJ   <- file.path(C01,  "cpm_BSJ_HnrnpM_PSD58.tsv")
CPM_FSJ   <- file.path(CFSJ, "cpm_FSJ_HnrnpM_PSD58.tsv")
OUT_DIR   <- file.path(ROOT, "JW/CircManuscript/Figures/Main")
OUT_PDF   <- file.path(OUT_DIR, "Fig_3C_v2.pdf")     # the ONLY file this script writes

# ---- Target physical size + uniform-shrink bookkeeping (matched to Fig_3C.R) ----
# Rendered at the SAME page size as Fig_3C.pdf (864 x 766 pt) so that dragging this PDF
# into PowerPoint and resizing it to 16 x 14.2 cm lands every inflated font at 8 pt,
# identical to Fig_3C.pdf. ggsave(scale = SC) enlarges the raw page by SC; fonts meant
# to end at 8 pt are pre-inflated to AXIS_PT * SC, and the gene-name strip keeps
# Fig_3C's original 6.6 * 1.5 pt proportion.
W_CM   <- 16; H_CM <- 14.2          # nominal FINAL size on the slide (same as Fig_3C.pdf)
ORIG_W_IN <- 12                     # Fig_3C rendered the grid at width = 12 in
SC     <- (ORIG_W_IN * 2.54) / W_CM # ggsave scale: draw big, shrink to target (= 1.905)
AXIS_PT <- 8                        # desired FINAL axis / title / legend / junction-strip font
axis_render_pt <- AXIS_PT * SC      # pre-inflate so it becomes 8 pt after the 1/SC shrink
MIN_PT <- 8                         # base_size (original proportion, as in Fig_3C.R)
NCOL   <- 9                         # circRNAs per block (36 = 4 x 9)

# ---- Sample order / palette (same convention as Fig_1E.R) ----
NEG_SAMPLES_ORD <- c("9_gNEG4_PSD5_S24","12_gNEG4_PSD5_S4","11_gNEG4_PSD8_S3","13_gNEG4_PSD8_S5")
HM_SAMPLES_ORD  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
SAMPLES_ORD     <- c(NEG_SAMPLES_ORD, HM_SAMPLES_ORD)
PAL_GROUP       <- c("gCTRL" = "#4C78A8", "gHNRNPM" = "#E45756")

# ---- Load data ----
rd <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "",
                             check.names = FALSE, stringsAsFactors = FALSE)
cls     <- rd(CLASS_ADD)
bsj_cpm <- rd(CPM_BSJ)
fsj_cpm <- rd(CPM_FSJ)

# ---- Select the 36 independent-back-splicing circRNAs (additive), most significant first ----
indep <- cls[cls$class == "independent_backsplicing", ]
indep <- indep[order(indep$adjP_BSJ), ]
cat("independent-back-splicing circRNAs:", nrow(indep), "\n")

# Unique gene-name panel labels (disambiguate duplicate host genes with #index)
gl  <- indep$gene_symbol
dup <- ave(seq_along(gl), gl, FUN = length) > 1
idx <- ave(seq_along(gl), gl, FUN = seq_along)
indep$panel_title <- ifelse(dup, paste0(gl, "#", idx), gl)

# Split the 36 into three stacked blocks of 12 (block 1 = most significant twelve)
indep$block <- ceiling(seq_len(nrow(indep)) / NCOL)

# ---- Long-format BSJ + FSJ CPM for those 36, carrying sample -> group + PSD ----
circs <- indep$circRNA
mk_long <- function(cpm_df, junction) {
  d <- cpm_df[match(circs, cpm_df$circRNA), ]        # order rows to match `circs`
  do.call(rbind, lapply(SAMPLES_ORD, function(s)
    data.frame(circRNA = circs, cpm = d[[s]], sample = s,
               junction = junction, stringsAsFactors = FALSE)))
}
pdat <- rbind(mk_long(bsj_cpm, "BSJ"), mk_long(fsj_cpm, "FSJ"))
pdat <- pdat[!is.na(pdat$cpm), ]

pdat$group    <- factor(ifelse(grepl("gHnrnpM", pdat$sample), "gHNRNPM", "gCTRL"),
                        levels = c("gCTRL","gHNRNPM"))
pdat$PSD      <- factor(sub(".*PSD([0-9]+).*", "\\1", pdat$sample), levels = c("5","8"))
pdat$junction <- factor(pdat$junction, levels = c("BSJ","FSJ"))
pdat$logCPM   <- log2(pdat$cpm + 1)
pdat$panel_title <- setNames(indep$panel_title, indep$circRNA)[pdat$circRNA]
pdat$block       <- setNames(indep$block,       indep$circRNA)[pdat$circRNA]

# Exactly three y-axis ticks per panel: placed at 10/50/90% of the panel's own data
# range (kept inside the range so rounding never drops one), rounded to 2 sig figs.
y3_breaks <- function(limits) {
  r <- range(limits, na.rm = TRUE)
  signif(r[1] + c(0.1, 0.5, 0.9) * diff(r), 2)
}

# ---- One circRNA = a single-column BSJ/FSJ mini-plot with its OWN independent y-axes ----
#   y-axis title + ticks are on the LEFT (ticks on every panel = independent y; title only
#   on the leftmost circRNA of a block, show_left). The BSJ/FSJ grey strips are on the
#   RIGHT and kept only on the rightmost circRNA of a block (show_right); the others drop
#   them. Fonts are pre-inflated (axis_render_pt) so they land at 8 pt after the 1/SC shrink.
make_circ <- function(circ, show_left, show_right) {
  sub <- pdat[pdat$circRNA == circ, ]
  g <- ggplot(sub, aes(x = group, y = logCPM, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey25", linewidth = 0.3) +
    geom_jitter(aes(shape = PSD), width = 0.16, size = 1.8, alpha = 0.85, colour = "grey20") +  # 2x the previous 0.9
    scale_fill_manual(values = PAL_GROUP, name = "Group") +
    scale_shape_manual(values = c(`5` = 16, `8` = 17), name = "PSD") +
    scale_y_continuous(breaks = y3_breaks,
                       labels = function(x) sprintf("%.2f", x)) +   # three ticks, all at 2 decimal places
    facet_grid(junction ~ panel_title, scales = "free_y") +         # junction strips on the RIGHT (default)
    labs(x = NULL, y = if (show_left) "log2(CPM + 1)" else NULL) +  # y title + ticks on the LEFT
    theme_bw(base_size = MIN_PT) +
    theme(
      text               = element_text(size = MIN_PT),
      axis.title         = element_text(size = axis_render_pt),      # -> 8 pt after shrink
      axis.text.x        = element_blank(),                          # group is redundant with the fill legend
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = axis_render_pt / 2),  # -> 4 pt after shrink (50% of original)
      strip.clip         = "off",                                     # let narrow-panel gene names overflow, not clip
      strip.text.x       = element_text(size = 6.6 * 1.5, face = "bold"),  # gene labels: Fig_3C proportion
      strip.text.y       = element_text(size = axis_render_pt, face = "bold", angle = 0),  # BSJ/FSJ -> 8 pt
      strip.background.x = element_rect(fill = "grey92", colour = "grey60"),
      strip.background.y = element_rect(fill = "grey80", colour = "grey60"),
      legend.position    = "bottom",
      legend.title       = element_text(size = axis_render_pt),
      legend.text        = element_text(size = axis_render_pt),
      panel.spacing.y    = grid::unit(0.30, "lines"))
  if (!show_right)
    g <- g + theme(strip.text.y       = element_blank(),
                   strip.background.y = element_blank())
  g
}

# ---- Assemble: 9 mini-plots per block (side by side), 4 blocks stacked, 1 bottom legend ----
build_block <- function(bl) {
  cs <- indep$circRNA[indep$block == bl]
  plots <- lapply(seq_along(cs), function(i)
    make_circ(cs[i], show_left = (i == 1), show_right = (i == length(cs))))
  wrap_plots(plots, nrow = 1)
}
block_rows <- lapply(sort(unique(indep$block)), build_block)

p <- wrap_plots(block_rows, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(OUT_PDF, p, width = W_CM, height = H_CM, units = "cm", scale = SC, device = cairo_pdf)
cat("Wrote:", OUT_PDF, "\n")
cat(sprintf("blocks: %d, columns/block: %d, circRNAs: %d (each panel independent y-axis)\n",
            length(block_rows), NCOL, nrow(indep)))

# ============================================================
# FIGURE LEGEND (for Fig_3C_v2.pdf)
# ============================================================
# Per-sample back-splice junction (BSJ) and forward-splice junction (FSJ) expression
# for the 36 independent-back-splicing circRNAs, the BSJ-DE circRNAs whose back-splice
# product changes under HNRNPM knockdown without a matching change in the linear (FSJ)
# junction (additive ~ group + PSD model; FSJ classed as unchanged under the soft rule
# of raw P >= 0.05 or |FSJ log2 fold change| < log2(1.5)). The 36 circRNAs are ordered
# by BSJ significance and arranged in four stacked blocks of nine columns; each
# column is one circRNA, labelled by its host gene symbol (duplicate host genes
# disambiguated with a #index). Within every column the upper sub-row (BSJ) and lower
# sub-row (FSJ) show log2(CPM + 1) for the control (gCTRL, left box) and HNRNPM
# knockdown (gHNRNPM, right box) groups, with individual samples overlaid as points
# whose shape encodes the post-selection day (PSD) and whose fill colour encodes the
# gRNA treatment group; the two groups are given by fill colour and the legend and are
# not re-labelled on the x-axis. Each panel is scaled on its own independent y-axis
# (tick labels on the left); the BSJ/FSJ strips sit at the right of each row. Boxes
# show the median and
# interquartile range. In this class the BSJ boxes shift between control and knockdown
# while the FSJ boxes at the same junction stay comparable, the signature of a change
# in back-splicing that is independent of the host gene's linear splicing.
