# Fig S3C v1 — the 20 NON-independent BSJ-DE circRNAs (co-regulated + opposite-
#              regulation + extreme-linear-collapse), drawn in the SAME style as
#              Fig_3C_v2.R. This is to Fig_S3C.R exactly what Fig_3C_v2.R is to Fig_3C.R:
#              same circRNA set, re-rendered in the Fig_1E-derived layout.
#
# Exports ONE file: Figures/Supplementary/Fig_S3C_v1.pdf
#
# KEPT FROM Fig_S3C.R
#   * the circRNA set: class in {co_regulated, opposite_regulation,
#     extreme_linear_collapse} under the additive (~ group + PSD) soft-rule
#     classification, ordered by class then adjP_BSJ (14 + 2 + 4 = 20);
#   * panel labels = host gene (#index for duplicates) PLUS the class abbreviation
#     (co-reg / opp / collapse), since the classes are mixed here;
#   * the same BSJ + FSJ CPM inputs and the same 7 samples (4 gCTRL + 3 gHNRNPM).
#
# ADOPTED FROM Fig_3C_v2.R (the new style)
#   * each circRNA is its OWN single-column facet_grid(junction ~ gene, scales="free_y")
#     mini-plot (BSJ sub-row over FSJ sub-row) => a FULLY INDEPENDENT y-axis per panel;
#     the mini-plots are laid side by side per block and stacked with patchwork;
#   * x = group (gCTRL vs gHNRNPM), 2-colour group fill, points shaped by PSD;
#   * y-axis tick labels on the LEFT (title once per block), exactly three ticks per
#     panel at 2 decimal places and half size; BSJ/FSJ grey strips on the RIGHT of each
#     row; x-axis group tick labels removed; single bottom legend; no title.
#
# LAYOUT: EIGHT circRNAs per row, stacked in three blocks (20 = 8 + 8 + 4); the short
#   last block is left-aligned with blank spacers so every panel keeps the same width.
#   Width is 16 cm with scale = SC and fonts pre-inflated to 8 pt (as Fig_3C_v2), so
#   resizing this PDF to 16 cm wide on the slide lands the fonts at 8 pt; the height
#   scales with the number of blocks. Note the 4 extreme-linear-collapse circRNAs have
#   no FSJ CPM, so their panels (all in the last block) show only a BSJ sub-row.

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
OUT_DIR   <- file.path(ROOT, "JW/CircManuscript/Figures/Supplementary")
OUT_PDF   <- file.path(OUT_DIR, "Fig_S3C_v1.pdf")    # the ONLY file this script writes

# ---- Target physical size + uniform-shrink bookkeeping (matched to Fig_3C_v2 / Fig_3C) ----
W_CM   <- 16                        # nominal FINAL width on the slide (same as Fig_3C_v2)
ROW_PER_BLOCK <- 14.2 / 4           # per-block height matched to Fig_3C_v2 (= 3.55 cm); H_CM set below
ORIG_W_IN <- 12                     # Fig_3C rendered the grid at width = 12 in
SC     <- (ORIG_W_IN * 2.54) / W_CM # ggsave scale: draw big, shrink to target (= 1.905)
AXIS_PT <- 8                        # desired FINAL axis / title / legend / junction-strip font
axis_render_pt <- AXIS_PT * SC      # pre-inflate so it becomes 8 pt after the 1/SC shrink
MIN_PT <- 8                         # base_size (original proportion, as in Fig_3C.R)
NCOL   <- 8                         # circRNAs per row (20 -> blocks of 8 + 8 + 4)

# ---- Sample order / palette (same convention as Fig_1E.R / Fig_3C_v2.R) ----
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

# ---- Select the 20 NON-independent BSJ-DE circRNAs (additive), grouped by class ----
other_classes <- c("co_regulated","opposite_regulation","extreme_linear_collapse")
other <- cls[cls$class %in% other_classes, ]
other <- other[order(match(other$class, other_classes), other$adjP_BSJ), ]   # class, then significance
cat("non-independent BSJ-DE circRNAs:", nrow(other), "\n")

# Panel labels: gene (#index for duplicate host genes) + class abbreviation
abbr <- c(co_regulated = "co-reg", opposite_regulation = "opp",
          extreme_linear_collapse = "collapse")
gl  <- other$gene_symbol
dup <- ave(seq_along(gl), gl, FUN = length) > 1
idx <- ave(seq_along(gl), gl, FUN = seq_along)
gene_label <- ifelse(dup, paste0(gl, "#", idx), gl)
other$panel_title <- sprintf("%s\n%s", gene_label, abbr[other$class])

# Split the 20 into stacked blocks of NCOL (8 + 8 + 4), in the class/significance order above
other$block <- ceiling(seq_len(nrow(other)) / NCOL)

# ---- Long-format BSJ + FSJ CPM for the 20, carrying sample -> group + PSD ----
circs <- other$circRNA
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
pdat$panel_title <- setNames(other$panel_title, other$circRNA)[pdat$circRNA]
pdat$block       <- setNames(other$block,       other$circRNA)[pdat$circRNA]

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
    geom_jitter(aes(shape = PSD), width = 0.16, size = 1.8, alpha = 0.85, colour = "grey20") +  # 2/3 of the previous 2.7
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
      strip.clip         = "off",                                     # let strip labels overflow, not clip
      strip.text.x       = element_text(size = 6.6 * 1.5, face = "bold"),  # gene+class labels: Fig_3C proportion
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

# ---- Assemble: up to NCOL mini-plots per block (side by side), blocks stacked, 1 bottom legend ----
# A short final block is padded with blank spacers so every panel keeps the same 1/NCOL
# width and the row is left-aligned (rather than stretched across the full width).
build_block <- function(bl) {
  cs <- other$circRNA[other$block == bl]
  plots <- lapply(seq_along(cs), function(i)
    make_circ(cs[i], show_left = (i == 1), show_right = (i == length(cs))))
  if (length(plots) < NCOL)
    plots <- c(plots, replicate(NCOL - length(plots), plot_spacer(), simplify = FALSE))
  wrap_plots(plots, nrow = 1)
}
block_rows <- lapply(sort(unique(other$block)), build_block)

H_CM <- ROW_PER_BLOCK * length(block_rows)     # height scales with the number of stacked blocks

p <- wrap_plots(block_rows, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(OUT_PDF, p, width = W_CM, height = H_CM, units = "cm", scale = SC, device = cairo_pdf)
cat("Wrote:", OUT_PDF, "\n")
cat(sprintf("blocks: %d, columns/block: %d, circRNAs: %d (each panel independent y-axis)\n",
            length(block_rows), NCOL, nrow(other)))

# ============================================================
# FIGURE LEGEND (for Fig_S3C_v1.pdf)
# ============================================================
# Per-sample back-splice junction (BSJ) and forward-splice junction (FSJ) expression
# for the 20 BSJ-DE circRNAs that are NOT independent back-splicing, i.e. the
# co-regulated, opposite-regulation, and extreme-linear-collapse classes under the
# additive ~ group + PSD model (soft FSJ rule). This is the companion to Fig_3C_v1/v2,
# which show the 36 independent circRNAs, rendered in the same style. The 20 circRNAs
# are ordered by class (co-regulated, then opposite-regulation, then
# extreme-linear-collapse) and within class by BSJ significance, and arranged eight per
# row in three stacked blocks (8 + 8 + 4); each column is one circRNA, labelled by its host gene
# symbol (duplicate host genes disambiguated with a #index) and its class abbreviation
# (co-reg / opp / collapse). Within every column the upper sub-row (BSJ) and lower
# sub-row (FSJ) show log2(CPM + 1) for the control (gCTRL, left box) and HNRNPM
# knockdown (gHNRNPM, right box) groups, with individual samples overlaid as points
# whose shape encodes the post-selection day (PSD) and whose fill colour encodes the
# gRNA treatment group; the two groups are given by fill colour and the legend and are
# not re-labelled on the x-axis. Each panel is scaled on its own independent y-axis
# (tick labels on the left); the BSJ/FSJ strips sit at the right of each row. Boxes
# show the median and interquartile range. Unlike the independent class, here the FSJ
# signal also changes (co-regulated / opposite) or is essentially absent (collapse),
# indicating the back-splicing change is not independent of the linear junction.
