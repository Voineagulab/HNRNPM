suppressPackageStartupMessages({
  library(SummarizedExperiment); library(limma); library(edgeR)
  library(ggplot2); library(patchwork)
})

# FinalAnalysis item 02 (gene): filterByExpr + TMM + voom + lmFit + eBayes.
# Two designs in parallel:
#   primary  : ~ group         (publication-primary)
#   additive : ~ group + PSD   (timepoint sensitivity comparator)
# No interaction model (n=1 HnrnpM_PSD8 -> uninterpretable).

GENE_RDS <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_rnaseq_230423/results/star_salmon/salmon.merged.gene_counts.rds"
OUT_DIR  <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/02_GeneTx_DE"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

HM  <- c("X5_gHnrnpM_PSD5_S20","X6_gHnrnpM_PSD5_S21","X7_gHnrnpM_PSD8_S22")
NEG <- c("X9_gNEG4_PSD5_S24","X11_gNEG4_PSD8_S3","X12_gNEG4_PSD5_S4","X13_gNEG4_PSD8_S5")
samps <- c(NEG, HM)

se <- readRDS(GENE_RDS)
counts_mat <- round(assay(se, "counts"))
annot      <- as.data.frame(rowData(se))
mat <- counts_mat[, samps]

group <- factor(c(rep("NEG4", length(NEG)), rep("HnrnpM", length(HM))),
                levels = c("NEG4","HnrnpM"))
PSD   <- factor(sub(".*PSD([0-9]+).*", "\\1", samps), levels = c("5","8"))
cat("Group:", as.character(group), "\nPSD:  ", as.character(PSD), "\n")

dge <- DGEList(counts = mat, group = group)
keep <- filterByExpr(dge, group = group)
cat("filterByExpr keeps", sum(keep), "/", length(keep),
    "genes (", round(100*sum(keep)/length(keep), 1), "%)\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

fit_de <- function(design, coef_name, label) {
  v   <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  out <- data.frame(gene_id = rownames(res),
                    gene_name = annot[rownames(res), "gene_name"],
                    res, check.names = FALSE)
  n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  n_up  <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_dn  <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
  cat(sprintf("[%s] coef '%s': tested %d, sig %d (UP %d, DOWN %d)\n",
              label, coef_name, nrow(out), n_sig, n_up, n_dn))
  list(res = out, v = v, fit = fit)
}

# Primary: ~ group
design_pri <- model.matrix(~ group); rownames(design_pri) <- samps
res_pri <- fit_de(design_pri, "groupHnrnpM", "PRIMARY ~group")
write.table(res_pri$res, file.path(OUT_DIR, "limma_gene_HnrnpM_vs_NEG_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Additive: ~ group + PSD
design_add <- model.matrix(~ group + PSD); rownames(design_add) <- samps
res_add <- fit_de(design_add, "groupHnrnpM", "ADDITIVE ~group+PSD")
write.table(res_add$res, file.path(OUT_DIR, "limma_gene_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
# PSD coefficient from the additive model
res_PSD <- topTable(res_add$fit, coef = "PSD8", number = Inf, sort.by = "P")
res_PSD_out <- data.frame(gene_id = rownames(res_PSD),
                          gene_name = annot[rownames(res_PSD), "gene_name"],
                          res_PSD, check.names = FALSE)
write.table(res_PSD_out, file.path(OUT_DIR, "limma_gene_PSDcoef_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
n_psd <- sum(res_PSD$adj.P.Val < 0.05, na.rm = TRUE)
cat(sprintf("[PSD coef] PSD8 vs PSD5: sig at adj.P<0.05 = %d\n", n_psd))

# ---- Comparison: primary vs additive ----
m <- merge(res_pri$res[, c("gene_id","logFC","adj.P.Val")],
           res_add$res[, c("gene_id","logFC","adj.P.Val")],
           by = "gene_id", suffixes = c("_primary","_add"))
m$sig_pri <- m$adj.P.Val_primary < 0.05
m$sig_add <- m$adj.P.Val_add     < 0.05
both <- sum(m$sig_pri & m$sig_add & sign(m$logFC_primary) == sign(m$logFC_add), na.rm = TRUE)
pri_only <- sum(m$sig_pri & !m$sig_add, na.rm = TRUE)
add_only <- sum(!m$sig_pri & m$sig_add, na.rm = TRUE)
flip <- sum(m$sig_pri & m$sig_add & sign(m$logFC_primary) != sign(m$logFC_add), na.rm = TRUE)
rho_all <- cor(m$logFC_primary, m$logFC_add, method = "spearman")
rho_sig <- cor(m$logFC_primary[m$sig_pri], m$logFC_add[m$sig_pri], method = "spearman")
cat(sprintf("\n--- Gene primary vs additive overlap ---\nprimary sig: %d\nadditive sig: %d\nboth (same dir): %d\nprimary-only: %d\nadditive-only: %d\nflipped: %d\nrho (all): %.3f\nrho (primary sig): %.3f\n",
            sum(m$sig_pri, na.rm=TRUE), sum(m$sig_add, na.rm=TRUE),
            both, pri_only, add_only, flip, rho_all, rho_sig))
write.table(data.frame(metric = c("primary_sig","additive_sig","both_sig_same_dir","primary_only","additive_only","flipped","spearman_logFC_all","spearman_logFC_primary_sig","universe","PSD_coef_sig"),
                       value  = c(sum(m$sig_pri,na.rm=TRUE), sum(m$sig_add,na.rm=TRUE), both, pri_only, add_only, flip, round(rho_all,4), round(rho_sig,4), nrow(res_pri$res), n_psd)),
            file.path(OUT_DIR, "gene_primary_vs_additive_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Volcano + scatter ----
volcano <- function(res, title_lab) {
  res$direction <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0, "UP in HnrnpM",
                   ifelse(res$adj.P.Val < 0.05 & res$logFC < 0, "DOWN in HnrnpM", "NS"))
  n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  ggplot(res, aes(x = logFC, y = -log10(P.Value), colour = direction)) +
    geom_point(alpha = 0.4, size = 0.6) +
    scale_colour_manual(values = c("UP in HnrnpM"="#d73027",
                                    "DOWN in HnrnpM"="#2166ac","NS"="grey75")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
    labs(title = title_lab,
         subtitle = sprintf("%d tested  ·  %d sig at adj.P<0.05", nrow(res), n_sig),
         x = "log2 Fold Change", y = "-log10(P)", colour = "") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}
ggsave(file.path(OUT_DIR, "volcano_gene_HnrnpM_PSD58.pdf"),
       volcano(res_pri$res, "limma-voom gene: ~group"), width = 8, height = 6)
ggsave(file.path(OUT_DIR, "volcano_gene_HnrnpM_PSD58_PSDadd.pdf"),
       volcano(res_add$res, "limma-voom gene: ~group + PSD"), width = 8, height = 6)

m$cat <- with(m, factor(ifelse(sig_pri & sig_add & sign(logFC_primary)==sign(logFC_add), "Both sig",
                       ifelse(sig_pri & !sig_add, "Primary only",
                       ifelse(!sig_pri & sig_add, "Additive only", "Neither/NS"))),
                       levels = c("Both sig","Primary only","Additive only","Neither/NS")))
pal <- c("Both sig"="#2ca02c","Primary only"="#d62728","Additive only"="#ff7f0e","Neither/NS"="grey80")
p_sc <- ggplot(m, aes(x = logFC_primary, y = logFC_add, colour = cat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey85") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey85") +
  geom_point(data = subset(m, cat=="Neither/NS"), size = 0.3, alpha = 0.2) +
  geom_point(data = subset(m, cat!="Neither/NS"), size = 0.6, alpha = 0.55) +
  scale_colour_manual(values = pal) +
  labs(title = "Gene logFC: primary vs additive (FinalAnalysis 02)",
       subtitle = sprintf("rho=%.3f overall, %.3f among primary sig", rho_all, rho_sig),
       x = "primary logFC (~group)", y = "additive logFC (~group+PSD)", colour = "") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "scatter_gene_primary_vs_additive.pdf"), p_sc, width = 8, height = 6.5)

cat("\nDone — gene DE (primary + additive).\n")
