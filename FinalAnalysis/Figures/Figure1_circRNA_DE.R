suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel)
  library(pheatmap); library(ggplotify)
  library(dplyr); library(tidyr)
})

# Figure 1: HnrnpM circRNA differential expression — publication figure
# All panels anchored on the ADDITIVE (~ group + PSD) model (publication primary).
#
#  A. Volcano (BSJ DE, additive)
#  B. CLIP enrichment per-region box (DE flank vs non-DE flank vs backgrounds)
#  C. CLIP enrichment per-circRNA mean flank score (DE vs non-DE)
#  D. Heatmap of 56 BSJ-DE x 7 samples (z-scored log2-CPM); row annotations = class + CLIP-bound
#  E. Top 10 independently-regulated BSJ-DE: per-hit boxplots (BSJ + FSJ x NEG + KO)

ROOT       <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis"
RESULTS    <- file.path(ROOT, "RESULTS")
OUT_DIR    <- file.path(ROOT, "Figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BSJ_ADD    <- file.path(RESULTS, "01_circRNA_DE", "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
CPM_BSJ    <- file.path(RESULTS, "01_circRNA_DE", "cpm_BSJ_HnrnpM_PSD58.tsv")
CPM_FSJ    <- file.path(RESULTS, "01_circRNA_DE", "FSJ", "cpm_FSJ_HnrnpM_PSD58.tsv")
CLASS_ADD  <- file.path(RESULTS, "01_circRNA_DE", "FSJ", "BSJ_FSJ_classification_additive.tsv")
CLIP_SCORE <- file.path(RESULTS, "09_HnrnpM_CLIP_enrichment", "circRNA_flank_score_PSD58.tsv")
CLIP_SIG   <- file.path(RESULTS, "09_HnrnpM_CLIP_enrichment", "clip_signal_per_region.tsv")
CLIP_TESTS <- file.path(RESULTS, "09_HnrnpM_CLIP_enrichment", "enrichment_test_summary.tsv")

NEG_SAMPLES_ORD <- c("9_gNEG4_PSD5_S24","12_gNEG4_PSD5_S4","11_gNEG4_PSD8_S3","13_gNEG4_PSD8_S5")
HM_SAMPLES_ORD  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
SAMPLES_ORD     <- c(NEG_SAMPLES_ORD, HM_SAMPLES_ORD)

# Palettes
PAL_DIR   <- c("UP in HnrnpM" = "#d73027", "DOWN in HnrnpM" = "#2166ac", "NS" = "grey75")
PAL_GROUP <- c("NEG4" = "#4C78A8", "HnrnpM" = "#E45756")
PAL_CLASS <- c("independent_backsplicing" = "#1f77b4",
               "co_regulated"             = "#2ca02c",
               "opposite_regulation"      = "#d62728",
               "extreme_linear_collapse"  = "#9467bd")

# ============================================================
# Load data
# ============================================================
bsj_add <- read.table(BSJ_ADD,    sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cls_add <- read.table(CLASS_ADD,  sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cpm_bsj <- read.table(CPM_BSJ,    sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cpm_fsj <- read.table(CPM_FSJ,    sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
clip    <- read.table(CLIP_SCORE, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
clip_sig   <- read.table(CLIP_SIG,   sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
clip_tests <- read.table(CLIP_TESTS, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

bsj_add$direction <- with(bsj_add,
  ifelse(adj.P.Val < 0.05 & logFC > 0, "UP in HnrnpM",
  ifelse(adj.P.Val < 0.05 & logFC < 0, "DOWN in HnrnpM", "NS")))
sig_set <- bsj_add[bsj_add$adj.P.Val < 0.05, ]
sig_ids <- sig_set$circRNA

# Merge classification (gives mechanistic class per sig circRNA)
sig_cls <- merge(sig_set[, c("circRNA","gene_symbol","logFC","adj.P.Val")],
                 cls_add[, c("circRNA","logFC_FSJ","adjP_FSJ","class")],
                 by="circRNA", all.x=TRUE)

cat("Loaded:\n",
    "  BSJ DE additive:   ", nrow(bsj_add), "(sig:", length(sig_ids), ")\n",
    "  Class additive:    ", nrow(cls_add), "\n",
    "  CLIP signal rows:  ", nrow(clip_sig), "\n",
    "  CLIP flank scores: ", nrow(clip), "\n")

# ============================================================
# PANEL A — Volcano
# ============================================================
n_sig <- sum(bsj_add$direction != "NS")
n_up  <- sum(bsj_add$direction == "UP in HnrnpM")
n_dn  <- sum(bsj_add$direction == "DOWN in HnrnpM")
top_lbl <- bsj_add[bsj_add$direction != "NS", ]
top_lbl <- top_lbl[order(top_lbl$adj.P.Val), ][1:10, ]

panel_A <- ggplot(bsj_add, aes(x = logFC, y = -log10(P.Value), colour = direction)) +
  geom_point(alpha = 0.55, size = 0.7) +
  scale_colour_manual(values = PAL_DIR, name = "") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey55") +
  geom_text_repel(data = top_lbl, aes(label = gene_symbol),
                  size = 2.6, min.segment.length = 0, max.overlaps = 20,
                  show.legend = FALSE, seed = 1) +
  labs(title = "A. Volcano: BSJ DE, additive (~ group + PSD)",
       subtitle = sprintf("%d tested · %d sig (UP %d / DOWN %d) at adj.P < 0.05",
                          nrow(bsj_add), n_sig, n_up, n_dn),
       x = "log2 Fold Change (HnrnpM vs NEG4)", y = "-log10(P-value)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey35"),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "Figure1A_volcano.pdf"), panel_A, width = 6, height = 5.5)

# ============================================================
# PANEL B — CLIP per-region log2(IP/Input) box plot (DE vs non-DE flank, backgrounds)
# ============================================================
de_circs <- sig_ids   # additive BSJ-DE
clip_sig$de_lab <- ifelse(clip_sig$circRNA %in% de_circs, "DE", "non-DE")
clip_sig$panel_label <- with(clip_sig, ifelse(
  region_type == "flank5", paste0(de_lab, " 5' flank"),
  ifelse(region_type == "flank3", paste0(de_lab, " 3' flank"),
  ifelse(region_type == "within_gene_bg", "Within-gene bg",
  ifelse(region_type == "non_host_bg",   "Non-host-gene bg", region_type)))))
clip_sig$panel_label <- factor(clip_sig$panel_label,
  levels = c("DE 5' flank","DE 3' flank","non-DE 5' flank","non-DE 3' flank",
             "Within-gene bg","Non-host-gene bg"))
clip_pal <- c("DE 5' flank" = "#d73027", "DE 3' flank" = "#d73027",
              "non-DE 5' flank" = "#4C78A8", "non-DE 3' flank" = "#4C78A8",
              "Within-gene bg" = "grey50", "Non-host-gene bg" = "grey75")

panel_B <- ggplot(clip_sig, aes(x = panel_label, y = log2_IP_over_Input, fill = panel_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, colour = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = clip_pal, guide = "none") +
  labs(title = "B. HnrnpM CLIP signal per region",
       subtitle = "DE = additive BSJ-DE (n=56); flanks of DE circRNAs vs flanks of non-DE circRNAs vs background introns",
       x = NULL, y = "log2(IP / Input)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

ggsave(file.path(OUT_DIR, "Figure1B_CLIP_per_region.pdf"), panel_B, width = 6, height = 5.5)

# ============================================================
# PANEL C — CLIP per-circRNA mean flank score (DE vs non-DE)
# ============================================================
de_n <- sum(clip$circRNA %in% de_circs)
nde_n <- sum(!(clip$circRNA %in% de_circs))
clip$is_DE_lab <- ifelse(clip$circRNA %in% de_circs,
                         sprintf("DE (n=%d)", de_n),
                         sprintf("non-DE (n=%d)", nde_n))
p_de_lab  <- sprintf("DE (n=%d)", de_n)
p_nde_lab <- sprintf("non-DE (n=%d)", nde_n)
p_anno <- clip_tests$p_one_sided[clip_tests$DE_axis == "BSJ_additive" &
                                   clip_tests$test == "1. DE flank vs non-DE flank"]

panel_C <- ggplot(clip, aes(x = is_DE_lab, y = mean_flank_log2, fill = is_DE_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_violin(alpha = 0.5, colour = "grey20", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  scale_fill_manual(values = setNames(c("#d73027","#4C78A8"), c(p_de_lab, p_nde_lab)),
                    guide = "none") +
  labs(title = "C. Per-circRNA mean flank score",
       subtitle = sprintf("Wilcoxon rank-sum (one-sided): P = %.2e", p_anno),
       x = NULL, y = "Mean log2(IP/Input) at flank introns") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"))

ggsave(file.path(OUT_DIR, "Figure1C_CLIP_per_circRNA.pdf"), panel_C, width = 5, height = 5.5)

# ============================================================
# PANEL D — Heatmap (56 BSJ-DE x 7 samples), row anns = Class + CLIP_bound
# ============================================================
cpm_mat <- as.matrix(cpm_bsj[match(sig_ids, cpm_bsj$circRNA), SAMPLES_ORD])
rownames(cpm_mat) <- sig_set$gene_symbol[match(sig_ids, sig_set$circRNA)]
dups <- duplicated(rownames(cpm_mat)) | duplicated(rownames(cpm_mat), fromLast = TRUE)
if (any(dups)) {
  rn <- rownames(cpm_mat)
  rn[dups] <- paste0(rn[dups], " (", sig_set$circRNA[match(sig_ids, sig_set$circRNA)][dups], ")")
  rownames(cpm_mat) <- rn
}

log_mat <- log2(cpm_mat + 1)
z_mat   <- t(scale(t(log_mat)))

col_ann <- data.frame(
  Group = factor(ifelse(grepl("gHnrnpM", SAMPLES_ORD), "HnrnpM", "NEG4"),
                 levels = c("NEG4","HnrnpM")),
  PSD = factor(sub(".*PSD([0-9]+).*", "\\1", SAMPLES_ORD), levels = c("5","8")),
  row.names = SAMPLES_ORD)

# Row annotations: only Class and CLIP_bound (per user request)
class_row <- cls_add$class[match(sig_ids, cls_add$circRNA)]
class_row[is.na(class_row)] <- "independent_backsplicing"

clip_row <- clip$bound_any_pos[match(sig_ids, clip$circRNA)]
clip_row[is.na(clip_row)] <- FALSE
clip_lab <- factor(ifelse(clip_row, "bound", "not bound"), levels = c("not bound","bound"))

row_ann <- data.frame(
  `CLIP_bound` = clip_lab,
  `Class`      = factor(class_row, levels = names(PAL_CLASS)),
  row.names = rownames(cpm_mat))

ann_colors <- list(
  Group      = c(NEG4 = "#4C78A8", HnrnpM = "#E45756"),
  PSD        = c(`5` = "#bdbdbd", `8` = "#525252"),
  CLIP_bound = c(bound = "#1a9850", `not bound` = "grey85"),
  Class      = PAL_CLASS
)

ph <- pheatmap(z_mat,
               cluster_rows = TRUE, cluster_cols = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_method = "average",
               annotation_col = col_ann,
               annotation_row = row_ann,
               annotation_colors = ann_colors,
               gaps_col = length(NEG_SAMPLES_ORD),
               color = colorRampPalette(c("#2166ac","white","#b2182b"))(101),
               breaks = seq(-2, 2, length.out = 102),
               border_color = NA,
               fontsize = 8, fontsize_row = 5.5, fontsize_col = 7,
               main = "D. 56 BSJ-DE circRNAs (additive) — log2(CPM+1) z-score per row",
               silent = TRUE)

ggsave(file.path(OUT_DIR, "Figure1D_heatmap.pdf"),
       as.ggplot(ph), width = 9, height = 9.5)

# ============================================================
# PANEL E — Top 10 INDEPENDENTLY-REGULATED circRNAs: per-hit boxplots
# ============================================================
indep <- sig_cls[sig_cls$class == "independent_backsplicing", ]
top10 <- indep[order(indep$adj.P.Val), ][1:10, ]
top10_ids <- top10$circRNA

cat("\nTop 10 independent_backsplicing circRNAs:\n")
print(top10[, c("circRNA","gene_symbol","logFC","adj.P.Val","logFC_FSJ","adjP_FSJ")], row.names=FALSE)

# Build long-format CPM tables
bsj_sub <- cpm_bsj[match(top10_ids, cpm_bsj$circRNA), c("circRNA","gene_symbol", SAMPLES_ORD)]
fsj_sub <- cpm_fsj[match(top10_ids, cpm_fsj$circRNA), c("circRNA","gene_symbol", SAMPLES_ORD)]

# Handle FSJ-filter-fail (FSJ missing → NA). Shouldn't happen for
# independent_backsplicing class but defensive.
missing_fsj <- top10_ids[is.na(fsj_sub$circRNA)]
if (length(missing_fsj) > 0) {
  filler <- data.frame(circRNA = missing_fsj,
                       gene_symbol = top10$gene_symbol[match(missing_fsj, top10_ids)],
                       check.names = FALSE)
  for (s in SAMPLES_ORD) filler[[s]] <- NA_real_
  fsj_sub[is.na(fsj_sub$circRNA), ] <- filler[order(match(filler$circRNA, top10_ids)), ]
}

bsj_long <- bsj_sub %>% pivot_longer(-c(circRNA, gene_symbol),
                                      names_to = "sample", values_to = "cpm") %>%
            mutate(junction = "BSJ")
fsj_long <- fsj_sub %>% pivot_longer(-c(circRNA, gene_symbol),
                                      names_to = "sample", values_to = "cpm") %>%
            mutate(junction = "FSJ")
long <- bind_rows(bsj_long, fsj_long)
long$group <- factor(ifelse(grepl("gHnrnpM", long$sample), "HnrnpM", "NEG4"),
                     levels = c("NEG4","HnrnpM"))
long$PSD <- factor(sub(".*PSD([0-9]+).*", "\\1", long$sample), levels = c("5","8"))
long$junction <- factor(long$junction, levels = c("BSJ","FSJ"))
long$logCPM   <- log2(long$cpm + 1)
gene_order_top10 <- top10$gene_symbol
long$gene_symbol <- factor(long$gene_symbol, levels = gene_order_top10)

panel_E <- ggplot(long, aes(x = group, y = logCPM, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey25", linewidth = 0.3) +
  geom_jitter(aes(shape = PSD), width = 0.16, size = 1.5,
              alpha = 0.85, colour = "grey20") +
  scale_fill_manual(values = PAL_GROUP, name = "Group") +
  scale_shape_manual(values = c(`5` = 16, `8` = 17), name = "PSD") +
  facet_grid(junction ~ gene_symbol, scales = "free_y", switch = "y") +
  labs(title = "E. Top 10 independently-regulated BSJ-DE circRNAs (additive) — per-sample CPM",
       subtitle = "Top row: BSJ counts (back-splice product); bottom row: FSJ counts (linear flux at same junction). For independent class: BSJ shifts, FSJ flat.",
       x = NULL, y = "log2(CPM + 1)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"),
        strip.background.x = element_rect(fill = "grey92", colour = "grey60"),
        strip.background.y = element_rect(fill = "grey80", colour = "grey60"),
        strip.text.x = element_text(size = 7.5, face = "bold"),
        strip.text.y.left = element_text(size = 8.5, face = "bold", angle = 0),
        axis.text.x = element_text(size = 7, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 7),
        legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        panel.spacing.x = unit(0.3, "lines"),
        panel.spacing.y = unit(0.3, "lines"))

ggsave(file.path(OUT_DIR, "Figure1E_top10_independent_boxplots.pdf"), panel_E, width = 14, height = 6.5)

# ============================================================
# Combined figure (5 panels: A|B|C on top row, D middle, E bottom)
# ============================================================
top_row    <- (panel_A | panel_B | panel_C) + plot_layout(widths = c(1.1, 1.1, 0.85))
mid_row    <- as.ggplot(ph)
bottom_row <- panel_E

combined <- (top_row / mid_row / bottom_row) +
  plot_layout(heights = c(1, 1.5, 1.1)) +
  plot_annotation(
    title = "Figure 1. HnrnpM-KD circRNA DE (additive model) and HnrnpM eCLIP enrichment at circRNA flanking introns",
    theme = theme(plot.title = element_text(size = 13, face = "bold")))

ggsave(file.path(OUT_DIR, "Figure1_circRNA_DE.pdf"),
       combined, width = 16, height = 22)

cat("\nWrote:\n",
    " ", file.path(OUT_DIR, "Figure1A_volcano.pdf"), "\n",
    " ", file.path(OUT_DIR, "Figure1B_CLIP_per_region.pdf"), "\n",
    " ", file.path(OUT_DIR, "Figure1C_CLIP_per_circRNA.pdf"), "\n",
    " ", file.path(OUT_DIR, "Figure1D_heatmap.pdf"), "\n",
    " ", file.path(OUT_DIR, "Figure1E_top10_independent_boxplots.pdf"), "\n",
    " ", file.path(OUT_DIR, "Figure1_circRNA_DE.pdf"), "\n")
