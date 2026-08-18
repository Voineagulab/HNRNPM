# Fig 3C v1 — The SAME 36 independent-back-splicing circRNAs as Fig_3C.R, but drawn
#             in the STYLE of Fig_1E.pdf (facet_grid junction x gene, group on the
#             x-axis, PSD encoded by point shape) instead of the 6 x 6 four-box grid.
#
# Exports ONE file: Figures/Main/Fig_3C_v1.pdf
#
# WHAT IS KEPT FROM Fig_3C.R
#   * the circRNA set: class == "independent_backsplicing" under the additive
#     (~ group + PSD) soft-rule classification, ordered by adjP_BSJ (most significant
#     first), 36 circRNAs, duplicate host genes disambiguated with a #index.
#   * the same BSJ + FSJ CPM inputs and the same 7 samples (4 gCTRL + 3 gHNRNPM).
#
# WHAT IS TAKEN FROM Fig_1E.R (the target style)
#   * per circRNA, BSJ on the top sub-row and FSJ on the bottom sub-row
#     (facet_grid(junction ~ gene), switch = "y");
#   * x-axis = group (gCTRL vs gHNRNPM), fill = group (2 colours, not 4);
#   * jittered sample points shaped by PSD (5 = circle, 8 = triangle);
#   * y = log2(CPM + 1), free-y per junction row, grey strip backgrounds, 8 pt fonts.
#
# LAYOUT CHANGE (requested)
#   Fig_1E puts 10 circRNAs in one row and the legend on the right. Here the 36
#   circRNAs are split into THREE blocks of 12 columns, stacked vertically; each block
#   is its own facet_grid(junction ~ gene) mini-panel (a BSJ sub-row + an FSJ sub-row
#   headed by its 12 gene labels). The three blocks are combined with patchwork and
#   share a SINGLE legend placed at the BOTTOM. The PDF stays 16 cm wide.

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
OUT_PDF   <- file.path(OUT_DIR, "Fig_3C_v1.pdf")     # the ONLY file this script writes

# ---- Target physical size on the slide ----
W_CM   <- 16; H_CM <- 22           # width x height as it should appear in PowerPoint
MIN_PT <- 8                         # base font size (points)
NCOL   <- 12                        # circRNAs per block (36 = 3 x 12)

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

# ---- One Fig_1E-style block (12 circRNAs, BSJ sub-row + FSJ sub-row) ----
make_block <- function(bl) {
  sub <- pdat[pdat$block == bl, ]
  sub$panel_title <- factor(sub$panel_title,
                            levels = indep$panel_title[indep$block == bl])
  ggplot(sub, aes(x = group, y = logCPM, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey25", linewidth = 0.3) +
    geom_jitter(aes(shape = PSD), width = 0.16, size = 0.9, alpha = 0.85, colour = "grey20") +
    scale_fill_manual(values = PAL_GROUP, name = "Group") +
    scale_shape_manual(values = c(`5` = 16, `8` = 17), name = "PSD") +
    facet_grid(junction ~ panel_title, scales = "free_y", switch = "y") +
    labs(x = NULL, y = "log2(CPM + 1)") +
    theme_bw(base_size = MIN_PT) +
    theme(
      text               = element_text(size = MIN_PT),
      axis.title         = element_text(size = MIN_PT),
      axis.text.x        = element_text(size = 6.5, angle = 30, hjust = 1),
      axis.text.y        = element_text(size = 7),
      strip.text.x       = element_text(size = 6.0, face = "bold"),   # gene-name column labels
      strip.text.y.left  = element_text(size = MIN_PT, face = "bold", angle = 0),
      strip.background.x = element_rect(fill = "grey92", colour = "grey60"),
      strip.background.y = element_rect(fill = "grey80", colour = "grey60"),
      legend.position    = "bottom",
      legend.title       = element_text(size = MIN_PT),
      legend.text        = element_text(size = MIN_PT),
      panel.spacing.x    = grid::unit(0.25, "lines"),
      panel.spacing.y    = grid::unit(0.30, "lines"))
}

blocks <- lapply(sort(unique(indep$block)), make_block)

# Stack the three blocks and collect their identical legends into one at the bottom.
p <- wrap_plots(blocks, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = sprintf("Independent back-splicing circRNAs (additive ~ group + PSD): BSJ vs FSJ CPM  (n = %d)",
                    nrow(indep)),
    theme = theme(plot.title = element_text(size = MIN_PT, hjust = 0.5))) &
  theme(legend.position = "bottom")

ggsave(OUT_PDF, p, width = W_CM, height = H_CM, units = "cm", device = cairo_pdf)
cat("Wrote:", OUT_PDF, "\n")
cat(sprintf("blocks: %d, columns/block: %d, circRNAs: %d\n",
            length(blocks), NCOL, nrow(indep)))

# ============================================================
# FIGURE LEGEND (for Fig_3C_v1.pdf)
# ============================================================
# Per-sample back-splice junction (BSJ) and forward-splice junction (FSJ) expression
# for the 36 independent-back-splicing circRNAs, the BSJ-DE circRNAs whose back-splice
# product changes under HNRNPM knockdown without a matching change in the linear (FSJ)
# junction (additive ~ group + PSD model; FSJ classed as unchanged under the soft rule
# of raw P >= 0.05 or |FSJ log2 fold change| < log2(1.5)). The 36 circRNAs are ordered
# by BSJ significance and arranged in three stacked blocks of twelve columns; each
# column is one circRNA, labelled by its host gene symbol (duplicate host genes
# disambiguated with a #index). Within every column the upper sub-row (BSJ) and lower
# sub-row (FSJ) show log2(CPM + 1) for the control (gCTRL) and HNRNPM knockdown
# (gHNRNPM) groups, with individual samples overlaid as points whose shape encodes the
# post-selection day (PSD) and whose fill colour encodes the gRNA treatment group.
# Boxes show the median and interquartile range. In this class the BSJ boxes shift
# between control and knockdown while the FSJ boxes at the same junction stay
# comparable, the signature of a change in back-splicing that is independent of the
# host gene's linear splicing.
