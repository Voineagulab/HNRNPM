suppressPackageStartupMessages({
  library(SummarizedExperiment); library(edgeR); library(limma)
  library(ggplot2); library(patchwork); library(ggrepel); library(matrixStats)
})

# FinalAnalysis 05 — 3-panel PCA on the new (FinalAnalysis) data:
#   circRNA  : log2(CPM+1) from CIRIquant BSJ counts (universe 5,820)
#   gene     : voom log2-CPM from salmon gene counts after filterByExpr
#   CLR      : per-circRNA, per-sample CLR (5,820 circRNAs)
# Samples coloured by group (NEG4 vs HnrnpM), shaped by PSD (5 / 8).

ROOT <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS"
CPM_F    <- file.path(ROOT, "01_circRNA_DE", "cpm_BSJ_HnrnpM_PSD58.tsv")
CLR_F    <- file.path(ROOT, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
GENE_RDS <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_rnaseq_230423/results/star_salmon/salmon.merged.gene_counts.rds"
OUT_PDF  <- file.path(ROOT, "05_Summary_Plots", "pca_3panel_HnrnpM_PSD58.pdf")

HM  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
NEG <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
samples <- c(NEG, HM)

do_pca <- function(mat, label) {
  rv <- rowVars(mat); mat <- mat[!is.na(rv) & rv > 0, , drop = FALSE]
  pr <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  v <- 100 * pr$sdev^2 / sum(pr$sdev^2)
  data.frame(sample = colnames(mat), PC1 = pr$x[,1], PC2 = pr$x[,2],
             label = label, pc1 = round(v[1],1), pc2 = round(v[2],1),
             n_features = nrow(mat))
}

# circRNA CPM
cpm_dat <- read.table(CPM_F, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
meta_circ <- c("circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id")
cpm_mat <- as.matrix(cpm_dat[, samples])
rownames(cpm_mat) <- cpm_dat$circRNA
pca_circ <- do_pca(log2(cpm_mat + 1), "circRNA log2(CPM+1)")

# Gene voom log2-CPM
se <- readRDS(GENE_RDS); counts <- round(assay(se, "counts"))
samples_X <- paste0("X", samples); mat_g <- counts[, samples_X]; colnames(mat_g) <- samples
group_g <- factor(ifelse(grepl("gHnrnpM", samples), "HnrnpM", "NEG4"),
                  levels = c("NEG4","HnrnpM"))
dge <- DGEList(counts = mat_g, group = group_g)
dge <- dge[filterByExpr(dge, group = group_g), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v   <- voom(dge, model.matrix(~group_g), plot = FALSE)
pca_gene <- do_pca(v$E, "gene log2-CPM (voom)")

# CLR
clr_dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
clr_mat <- as.matrix(clr_dat[, samples]); rownames(clr_mat) <- clr_dat$circRNA
pca_clr <- do_pca(clr_mat, "CLR (5,820 circRNAs)")

build <- function(d) {
  d$group <- factor(ifelse(grepl("gHnrnpM", d$sample), "HnrnpM", "NEG4"),
                    levels = c("NEG4","HnrnpM"))
  d$PSD   <- factor(sub(".*PSD([0-9]+).*", "\\1", d$sample), levels = c("5","8"))
  d$short <- sub("^[0-9]+_g", "", d$sample)
  ggplot(d, aes(x = PC1, y = PC2, colour = group, shape = PSD, label = short)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(size = 2.6, min.segment.length = 0, seed = 1, show.legend = FALSE) +
    scale_colour_manual(values = c(NEG4 = "#4C78A8", HnrnpM = "#E45756")) +
    scale_shape_manual(values = c(`5` = 16, `8` = 17)) +
    labs(title = unique(d$label),
         subtitle = sprintf("%d features", unique(d$n_features)),
         x = sprintf("PC1 (%.1f%%)", unique(d$pc1)),
         y = sprintf("PC2 (%.1f%%)", unique(d$pc2))) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}
p1 <- build(pca_circ); p2 <- build(pca_gene); p3 <- build(pca_clr)
combined <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "FinalAnalysis 05 — PCA: NEG4 vs HnrnpM (PSD5+8 only, 7 samples)",
    subtitle = sprintf("circRNA log2(CPM+1) (n=%d) | gene log2-CPM voom (n=%d) | CLR (n=%d)",
                       nrow(cpm_mat), nrow(v$E), nrow(clr_mat)),
    theme = theme(plot.title = element_text(size = 12),
                  plot.subtitle = element_text(size = 10, colour = "grey40"),
                  legend.position = "bottom")
  )
ggsave(OUT_PDF, combined, width = 13, height = 5.8)
cat("Wrote:", OUT_PDF, "\n")
