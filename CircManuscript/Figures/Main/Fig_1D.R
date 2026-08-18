# Fig 1D — Heatmap of the 56 additive BSJ-DE circRNAs x 7 samples.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/Figures/Figure1_circRNA_DE.R
# Retains only PANEL D of that figure — pheatmap of z-scored log2(CPM+1),
#   rows = sig circRNAs (gene symbol), cols = samples (NEG4 | HnrnpM, gap between),
#   row annotations = mechanistic Class + CLIP-bound; col annotations = Group + PSD.
# DE = additive BSJ-DE circRNAs (~ group + PSD, adj.P<0.05, publication primary).
# Saved as PDF sized for the manuscript slide (W_CM x H_CM), fonts 8 pt, via
#   cairo_pdf, to Figures/Main/Fig_1D.pdf.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(pheatmap); library(ggplotify)
})


# --- Path definition ---
ROOT       <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES        <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS")
BSJ_ADD    <- file.path(RES, "01_circRNA_DE", "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
CPM_BSJ    <- file.path(RES, "01_circRNA_DE", "cpm_BSJ_HnrnpM_PSD58.tsv")
CLASS_ADD  <- file.path(RES, "01_circRNA_DE", "FSJ", "BSJ_FSJ_classification_additive.tsv")
CLIP_SCORE <- file.path(RES, "09_HnrnpM_CLIP_enrichment", "circRNA_flank_score_PSD58.tsv")
OUT_DIR    <- file.path(ROOT, "JW/CircManuscript/Figures/Main")


# ---- Sample order (NEG4 first, then HnrnpM; gap drawn between groups) ----
NEG_SAMPLES_ORD <- c("9_gNEG4_PSD5_S24","12_gNEG4_PSD5_S4","11_gNEG4_PSD8_S3","13_gNEG4_PSD8_S5")
HM_SAMPLES_ORD  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
SAMPLES_ORD     <- c(NEG_SAMPLES_ORD, HM_SAMPLES_ORD)

PAL_CLASS <- c("independent_backsplicing" = "#1f77b4",
               "co_regulated"             = "#2ca02c",
               "opposite_regulation"      = "#d62728",
               "extreme_linear_collapse"  = "#9467bd")

# ---- Load data ----
bsj_add <- read.table(BSJ_ADD,    sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cls_add <- read.table(CLASS_ADD,  sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cpm_bsj <- read.table(CPM_BSJ,    sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
clip    <- read.table(CLIP_SCORE, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

sig_set <- bsj_add[bsj_add$adj.P.Val < 0.05, ]
sig_ids <- sig_set$circRNA
cat("sig circRNAs:", length(sig_ids), "\n")

# ---- Build z-scored matrix (rows = gene symbols, cols = samples) ----
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

# ---- Annotations ----
col_ann <- data.frame(
  Group = factor(ifelse(grepl("gHnrnpM", SAMPLES_ORD), "gHNRNPM", "gCTRL"),
                 levels = c("gCTRL","gHNRNPM")),
  PSD = factor(sub(".*PSD([0-9]+).*", "\\1", SAMPLES_ORD), levels = c("5","8")),
  row.names = SAMPLES_ORD)

class_row <- cls_add$class[match(sig_ids, cls_add$circRNA)]
class_row[is.na(class_row)] <- "independent_backsplicing"

clip_row <- clip$bound_any_pos[match(sig_ids, clip$circRNA)]
clip_row[is.na(clip_row)] <- FALSE
clip_lab <- factor(ifelse(clip_row, "bound", "not bound"), levels = c("not bound","bound"))

# Short single-line Class labels. NOTE: pheatmap allots only ONE line of height
# per annotation-legend entry, so multi-line ("\n") labels overlap the entry
# below. Keep these single-line; shorten instead of wrapping.
class_disp <- c(independent_backsplicing = "Indep",
                co_regulated             = "Co-reg",
                opposite_regulation      = "Opposite",
                extreme_linear_collapse  = "Collapse")

row_ann <- data.frame(
  `eCLIP` = clip_lab,
  `Class`        = factor(class_disp[class_row], levels = class_disp[names(PAL_CLASS)]),
  row.names = rownames(cpm_mat),
  check.names = FALSE)               # keep the hyphen in "eCLIP"

ann_colors <- list(
  Group          = c(gCTRL = "#4C78A8", gHNRNPM = "#E45756"),
  PSD            = c(`5` = "#bdbdbd", `8` = "#525252"),
  `eCLIP` = c(bound = "#1a9850", `not bound` = "grey85"),
  Class          = setNames(unname(PAL_CLASS), class_disp[names(PAL_CLASS)])
)

# Column display labels: rename groups (gNEG4 -> gCTRL, gHnrnpM -> gHNRNPM).
# SAMPLES_ORD itself is left unchanged — it must match the data column names.
labels_col <- unname(c("9_gNEG4_PSD5_S24" = "gCTRL_PSD5_rep1", "12_gNEG4_PSD5_S4" = "gCTRL_PSD5_rep2", "11_gNEG4_PSD8_S3" = "gCTRL_PSD8_rep1", "13_gNEG4_PSD8_S5" = "gCTRL_PSD8_rep2", "5_gHnrnpM_PSD5_S20" = "gHNRNPM_PSD5_rep1", "6_gHnrnpM_PSD5_S21" = "gHNRNPM_PSD5_rep2", "7_gHnrnpM_PSD8_S22" = "gHNRNPM_PSD8")[SAMPLES_ORD])

# ---- Heatmap: with no labels ----
W_CM   <- 16; H_CM <- 9.55         # width x height as it should appear in PowerPoint
MIN_PT <- 4                        
ROW_PT <- 4                        
LEG_PT <- 4                       

ph_nolegend <- pheatmap(z_mat,
               cluster_rows = TRUE, cluster_cols = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_method = "average",
               annotation_col = col_ann,
               annotation_row = row_ann,
               annotation_colors = ann_colors,
               gaps_col = length(NEG_SAMPLES_ORD),
               labels_col = labels_col,
               angle_col = "0",
               color = colorRampPalette(c("#2166ac","white","#b2182b"))(101),
               breaks = seq(-2, 2, length.out = 102),
               border_color = NA,
               fontsize = LEG_PT, fontsize_col = MIN_PT,
               fontsize_row = ROW_PT, 
               main = NA,
               silent = TRUE,
               legend = FALSE,             # hides the colour-scale legend
               annotation_legend = FALSE,  # hides the annotation (Group/PSD/eCLIP/Class) legends
               show_colnames = FALSE,
               annotation_names_col = FALSE,
               annotation_names_row = FALSE
               )

ggsave(file.path(OUT_DIR, "Fig_1D_nolegend.pdf"), 
       as.ggplot(ph_nolegend) +
         theme(plot.margin = margin(t = 2, r = 2, b = 14, l = 2, unit = "pt")),
       width = W_CM, height = H_CM, units = "cm",
       device = cairo_pdf)



# ---- Heatmap: with labels ----
W_CM   <- 16; H_CM <- 9.55          # width x height as it should appear in PowerPoint
MIN_PT <- 10                        # font size (points): column (sample) labels
LEG_PT <- 10                       # Must be big enough to avoid overlapping

ph_legend <- pheatmap(z_mat,
                cluster_rows = TRUE, cluster_cols = FALSE,
                clustering_distance_rows = "euclidean",
                clustering_method = "average",
                annotation_col = col_ann,
                annotation_row = row_ann,
                annotation_colors = ann_colors,
                gaps_col = length(NEG_SAMPLES_ORD),
                labels_col = labels_col,
                angle_col = "0",
                color = colorRampPalette(c("#2166ac","white","#b2182b"))(101),
                breaks = seq(-2, 2, length.out = 102),
                border_color = NA,
                fontsize_row = ROW_PT,
                fontsize = LEG_PT, fontsize_col = MIN_PT,
                main = "56 BSJ-DE circRNAs (additive) - log2(CPM+1) z-score per row.pdf",
                silent = TRUE)
ggsave(file.path(OUT_DIR, "Fig_1D_legend.pdf"), as.ggplot(ph_legend) +
         theme(plot.margin = margin(t = 2, r = 2, b = 14, l = 2, unit = "pt")),
       width = W_CM, height = H_CM, units = "cm",
       device = cairo_pdf) # Zoom to 6.71cm (height) x 11.25cm (width)
                        


# ---- Heatmap: with labels and vertical ----
W_CM   <- 13; H_CM <- 14.53         # width x height as it should appear in PowerPoint
MIN_PT <- 6                        
ROW_PT <- 6                        
LEG_PT <- 6                       

ph_nolegend <- pheatmap(z_mat,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        clustering_distance_rows = "euclidean",
                        clustering_method = "average",
                        annotation_col = col_ann,
                        annotation_row = row_ann,
                        annotation_colors = ann_colors,
                        gaps_col = length(NEG_SAMPLES_ORD),
                        labels_col = labels_col,
                        angle_col = "90",
                        color = colorRampPalette(c("#2166ac","white","#b2182b"))(101),
                        breaks = seq(-2, 2, length.out = 102),
                        border_color = NA,
                        fontsize = LEG_PT, fontsize_col = MIN_PT,
                        fontsize_row = ROW_PT, 
                        main = NA,
                        silent = TRUE,
                        legend = FALSE,             # hides the colour-scale legend
                        annotation_legend = FALSE,  # hides the annotation (Group/PSD/eCLIP/Class) legends
                        # show_colnames = FALSE,
                        # annotation_names_col = FALSE,
                        # annotation_names_row = FALSE
)

ggsave(file.path(OUT_DIR, "Fig_1D_nolegend_vertical.pdf"), 
       as.ggplot(ph_nolegend) +
         theme(plot.margin = margin(t = 2, r = 2, b = 0, l = 2, unit = "pt")),
       width = W_CM, height = H_CM, units = "cm",
       device = cairo_pdf)
                        
# # ---- Figure legend ----
# (D) Heatmap of the 56 circRNAs significantly differentially expressed at the BSJ level under the additive model at adjusted P < 0.05, shown across the seven samples. Colour encodes the per row z score of log2 of BSJ counts per million plus one, with red indicating expression above that circRNA row mean and blue below the row mean, and white near the row mean, with values clamped to the range minus two to plus two. Rows are circRNAs labelled by host gene symbol and clustered by Euclidean distance with average linkage, and columns are the individual samples in fixed order with a gap separating the four control samples (gCTRL) from the three HNRNPM knockdown samples (gHNRNPM). The column annotation bars indicate experimental group (gCTRL or gHNRNPM) and postnatal day (5 or 8). The row annotation bars indicate the mechanistic class of each circRNA (Indep for independent backsplicing, Co-reg for co-regulated, Opposite for opposite regulation, and Collapse for extreme linear collapse) and whether the flanking introns are bound in HNRNPM eCLIP (bound or not bound).
