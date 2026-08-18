# Fig 1E — Top 10 independently-regulated BSJ-DE circRNAs: per-sample CPM boxplots.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/Figures/Figure1_circRNA_DE.R
# Retains only PANEL E of that figure — BSJ (top) + FSJ (bottom) log2(CPM+1)
#   boxplots for the 10 most significant "independent_backsplicing" circRNAs,
#   faceted junction x gene, split by group (NEG4 vs HnrnpM). Panels A/B/C/D and
#   the combined assembly are dropped.
# DE = additive BSJ-DE circRNAs (~ group + PSD, adj.P<0.05, publication primary).
# Saved as PDF sized for the manuscript slide (W_CM x H_CM), fonts 8 pt, via
#   cairo_pdf, to Figures/Main/Fig_1E.pdf.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})


# --- Path definition ---
ROOT      <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES       <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
BSJ_ADD   <- file.path(RES, "01_circRNA_DE", "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
CPM_BSJ   <- file.path(RES, "01_circRNA_DE", "cpm_BSJ_HnrnpM_PSD58.tsv")
CPM_FSJ   <- file.path(RES, "01_circRNA_DE", "FSJ", "cpm_FSJ_HnrnpM_PSD58.tsv")
CLASS_ADD <- file.path(RES, "01_circRNA_DE", "FSJ", "BSJ_FSJ_classification_additive.tsv")
OUT_DIR   <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size on the slide ----
W_CM   <- 16; H_CM <- 6.83         # width x height as it should appear in PowerPoint
MIN_PT <- 8                         # font size (points) — applied to all text

NEG_SAMPLES_ORD <- c("9_gNEG4_PSD5_S24","12_gNEG4_PSD5_S4","11_gNEG4_PSD8_S3","13_gNEG4_PSD8_S5")
HM_SAMPLES_ORD  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
SAMPLES_ORD     <- c(NEG_SAMPLES_ORD, HM_SAMPLES_ORD)
PAL_GROUP       <- c("gCTRL" = "#4C78A8", "gHNRNPM" = "#E45756")

# ---- Load data ----
bsj_add <- read.table(BSJ_ADD,   sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cls_add <- read.table(CLASS_ADD, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cpm_bsj <- read.table(CPM_BSJ,   sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cpm_fsj <- read.table(CPM_FSJ,   sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

sig_set <- bsj_add[bsj_add$adj.P.Val < 0.05, ]
sig_cls <- merge(sig_set[, c("circRNA","gene_symbol","logFC","adj.P.Val")],
                 cls_add[, c("circRNA","logFC_FSJ","adjP_FSJ","class")],
                 by="circRNA", all.x=TRUE)

# ---- Top 10 independently-regulated circRNAs by significance ----
indep <- sig_cls[sig_cls$class == "independent_backsplicing", ]
top10 <- indep[order(indep$adj.P.Val), ][1:10, ]
top10_ids <- top10$circRNA
cat("Top 10 independent_backsplicing circRNAs:\n")
print(top10[, c("circRNA","gene_symbol","logFC","adj.P.Val")], row.names=FALSE)

# ---- Long-format BSJ + FSJ CPM ----
bsj_sub <- cpm_bsj[match(top10_ids, cpm_bsj$circRNA), c("circRNA","gene_symbol", SAMPLES_ORD)]
fsj_sub <- cpm_fsj[match(top10_ids, cpm_fsj$circRNA), c("circRNA","gene_symbol", SAMPLES_ORD)]

# Defensive: fill any FSJ-filter-fail rows with NA (shouldn't occur for this class)
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
long$group     <- factor(ifelse(grepl("gHnrnpM", long$sample), "gHNRNPM", "gCTRL"),
                         levels = c("gCTRL","gHNRNPM"))
long$PSD       <- factor(sub(".*PSD([0-9]+).*", "\\1", long$sample), levels = c("5","8"))
long$junction  <- factor(long$junction, levels = c("BSJ","FSJ"))
long$logCPM    <- log2(long$cpm + 1)
long$gene_symbol <- factor(long$gene_symbol, levels = top10$gene_symbol)

# ---- Panel E ----
pe <- ggplot(long, aes(x = group, y = logCPM, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey25", linewidth = 0.3) +
  geom_jitter(aes(shape = PSD), width = 0.16, size = 1.0, alpha = 0.85, colour = "grey20") +
  scale_fill_manual(values = PAL_GROUP, name = "Group") +
  scale_shape_manual(values = c(`5` = 16, `8` = 17), name = "PSD") +
  facet_grid(junction ~ gene_symbol, scales = "free_y", switch = "y") +
  labs(x = NULL, y = "log2(CPM + 1)") +
  theme_bw(base_size = MIN_PT) +
  theme(
    text              = element_text(size = MIN_PT),
    axis.title        = element_text(size = MIN_PT),
    axis.text.x       = element_text(size = MIN_PT, angle = 30, hjust = 1),
    axis.text.y       = element_text(size = MIN_PT),
    strip.text.x      = element_text(size = 6.5, face = "bold"),   # gene-name column labels
    strip.text.y.left = element_text(size = MIN_PT, face = "bold", angle = 0),
    strip.background.x = element_rect(fill = "grey92", colour = "grey60"),
    strip.background.y = element_rect(fill = "grey80", colour = "grey60"),
    legend.position   = "right",
    legend.title      = element_text(size = MIN_PT),
    legend.text       = element_text(size = MIN_PT),
    panel.spacing.x   = grid::unit(0.3, "lines"),
    panel.spacing.y   = grid::unit(0.3, "lines"))

ggsave(file.path(OUT_DIR, "Fig_1E.pdf"), pe,
       width = W_CM, height = H_CM, units = "cm",
       device = cairo_pdf)

# # ---- Figure legend ----
# (E) Per-sample expression box plots for the ten most significant independently regulated circRNAs, defined as circRNAs in the independent backsplicing class ranked by adjusted P value, with back-splice junction (BSJ) counts shown on top and forward-splice junction (FSJ) on the bottom. The ten circRNAs in columns labelled by host gene symbol. Boxes show the median and interquartile range, and individual points show the seven samples, with point shape encoding PSDs and fill colour encoding gRNA treatment group. 