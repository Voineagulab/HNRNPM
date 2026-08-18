suppressPackageStartupMessages({
  library(limma); library(edgeR); library(ggplot2); library(patchwork)
})

# FinalAnalysis item 01: circRNA differential expression from CIRIquant BSJ counts.
# 7 samples (3 HnrnpM_PSD5/8 vs 4 NEG4_PSD5/8). Two models fitted in parallel:
#   primary  : ~ group         (publication-primary)
#   additive : ~ group + PSD   (timepoint sensitivity; comparator)
# Library-size normalisation uses MultiQC R1 reads (preferred for low-count
# circRNA data over DESeq2 median-of-ratios).

ROOT     <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
COUNTS_F <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/bsj_matrix_CIRIquant_HnrnpM_PSD58.tsv")
MULTIQC  <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_circrna_230420/results/multiqc/multiqc_data/multiqc_general_stats.txt"
OUT_DIR  <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE")

ANNOT_COLS <- c("gene_symbol","chr","start","end","strand","circ_type","gene_id")

raw <- read.table(COUNTS_F, header = TRUE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE)
annot   <- raw[, c("circRNA", ANNOT_COLS), drop = FALSE]
counts  <- round(as.matrix(raw[, !(colnames(raw) %in% c("circRNA", ANNOT_COLS))]))
rownames(counts) <- raw$circRNA
samples <- colnames(counts)
cat("Loaded counts:", nrow(counts), "circRNAs x", ncol(counts), "samples\n")

group <- factor(ifelse(grepl("gHnrnpM", samples), "gHnrnpM", "gNEG4"),
                levels = c("gNEG4","gHnrnpM"))
PSD   <- factor(sub(".*PSD([0-9]+).*", "\\1", samples), levels = c("5","8"))
cat("Group:", as.character(group), "\nPSD:  ", as.character(PSD), "\n")

# Library-size from MultiQC R1
stats <- read.table(MULTIQC, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
r1    <- stats[grepl("_1$", rownames(stats)),
               "FastQC_mqc-generalstats-fastqc-total_sequences", drop = FALSE]
rownames(r1) <- sub("_1$", "", rownames(r1))
lib_sizes <- r1[samples, 1]; names(lib_sizes) <- samples
cat("Library sizes (R1 reads):\n"); print(lib_sizes)

dge <- DGEList(counts = counts, lib.size = lib_sizes, group = group)

# CPM table
cpm_mat <- cpm(dge)
cpm_out <- cbind(circRNA = rownames(cpm_mat), annot[match(rownames(cpm_mat), annot$circRNA), ANNOT_COLS],
                 as.data.frame(cpm_mat))
write.table(cpm_out, file.path(OUT_DIR, "cpm_BSJ_HnrnpM_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Helper: fit a design, extract a coefficient, write TSV, return result df
fit_de <- function(design, coef_name, label) {
  v   <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  out <- cbind(circRNA = rownames(res),
               annot[match(rownames(res), annot$circRNA), ANNOT_COLS],
               res)
  n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  n_up  <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_dn  <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
  cat(sprintf("[%s] coef '%s': tested %d, sig %d (UP %d, DOWN %d)\n",
              label, coef_name, nrow(out), n_sig, n_up, n_dn))
  return(out)
}

# ---- Primary: ~ group ----
design_pri <- model.matrix(~ group); rownames(design_pri) <- samples
res_pri <- fit_de(design_pri, "groupgHnrnpM", "PRIMARY ~group")
write.table(res_pri, file.path(OUT_DIR, "limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Additive: ~ group + PSD ----
design_add <- model.matrix(~ group + PSD); rownames(design_add) <- samples
res_add <- fit_de(design_add, "groupgHnrnpM", "ADDITIVE ~group+PSD")
write.table(res_add, file.path(OUT_DIR, "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Also save the PSD8 coefficient table from the additive model
v_add   <- voom(dge, design_add, plot = FALSE)
fit_add <- eBayes(lmFit(v_add, design_add))
res_PSD <- topTable(fit_add, coef = "PSD8", number = Inf, sort.by = "P")
res_PSD_out <- cbind(circRNA = rownames(res_PSD),
                      annot[match(rownames(res_PSD), annot$circRNA), ANNOT_COLS],
                      res_PSD)
write.table(res_PSD_out, file.path(OUT_DIR, "limma_BSJ_PSDcoef_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
n_psd <- sum(res_PSD$adj.P.Val < 0.05, na.rm = TRUE)
cat(sprintf("[PSD coef] PSD8 vs PSD5: sig at adj.P<0.05 = %d\n", n_psd))

# ---- Comparison: primary vs additive ----
m <- merge(res_pri[, c("circRNA","logFC","adj.P.Val")],
           res_add[, c("circRNA","logFC","adj.P.Val")],
           by = "circRNA", suffixes = c("_primary","_add"))
m$sig_pri <- m$adj.P.Val_primary < 0.05
m$sig_add <- m$adj.P.Val_add     < 0.05
both <- sum(m$sig_pri & m$sig_add & (sign(m$logFC_primary) == sign(m$logFC_add)), na.rm = TRUE)
pri_only <- sum(m$sig_pri & !m$sig_add, na.rm = TRUE)
add_only <- sum(!m$sig_pri & m$sig_add, na.rm = TRUE)
flip <- sum(m$sig_pri & m$sig_add & (sign(m$logFC_primary) != sign(m$logFC_add)), na.rm = TRUE)
rho_all <- cor(m$logFC_primary, m$logFC_add, method = "spearman")
rho_sig <- if (sum(m$sig_pri, na.rm=TRUE) >= 3)
              cor(m$logFC_primary[m$sig_pri], m$logFC_add[m$sig_pri], method = "spearman") else NA_real_
cat(sprintf("\n--- Primary vs additive overlap ---\n"))
cat(sprintf("primary sig: %d\nadditive sig: %d\nboth (same dir): %d\nprimary-only: %d\nadditive-only: %d\nflipped: %d\nrho (all): %.3f\nrho (primary sig): %.3f\n",
            sum(m$sig_pri, na.rm=TRUE), sum(m$sig_add, na.rm=TRUE),
            both, pri_only, add_only, flip, rho_all, rho_sig))

ov <- data.frame(
  metric = c("primary_sig","additive_sig","both_sig_same_dir","primary_only","additive_only",
             "sig_both_flipped","spearman_logFC_all","spearman_logFC_primary_sig",
             "universe","PSD_coef_sig"),
  value  = c(sum(m$sig_pri, na.rm=TRUE), sum(m$sig_add, na.rm=TRUE),
             both, pri_only, add_only, flip,
             round(rho_all, 4), round(rho_sig, 4),
             nrow(res_pri), n_psd))
write.table(ov, file.path(OUT_DIR, "BSJ_primary_vs_additive_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Volcano plots ----
volcano <- function(res, title_lab) {
  res$direction <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0, "UP in HnrnpM",
                   ifelse(res$adj.P.Val < 0.05 & res$logFC < 0, "DOWN in HnrnpM", "NS"))
  n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  n_up  <- sum(res$direction == "UP in HnrnpM")
  n_dn  <- sum(res$direction == "DOWN in HnrnpM")
  ggplot(res, aes(x = logFC, y = -log10(P.Value), colour = direction)) +
    geom_point(alpha = 0.5, size = 0.9) +
    scale_colour_manual(values = c("UP in HnrnpM"="#d73027",
                                    "DOWN in HnrnpM"="#2166ac","NS"="grey75")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
    labs(title = title_lab,
         subtitle = sprintf("%d tested  ·  %d sig (UP %d / DOWN %d) at adj.P<0.05",
                            nrow(res), n_sig, n_up, n_dn),
         x = "log2 Fold Change", y = "-log10(P)", colour = "") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}
p_pri <- volcano(res_pri, "limma-voom BSJ: ~group  (FinalAnalysis 01)")
p_add <- volcano(res_add, "limma-voom BSJ: ~group + PSD  (FinalAnalysis 01)")
ggsave(file.path(OUT_DIR, "volcano_BSJ_HnrnpM_PSD58.pdf"),         p_pri, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "volcano_BSJ_HnrnpM_PSD58_PSDadd.pdf"), p_add, width = 8, height = 6)

# ---- Scatter: primary vs additive logFC ----
m$cat <- with(m, factor(ifelse(sig_pri & sig_add & sign(logFC_primary) == sign(logFC_add), "Both sig",
                       ifelse(sig_pri & !sig_add, "Primary only",
                       ifelse(!sig_pri & sig_add, "Additive only", "Neither/NS"))),
                       levels = c("Both sig","Primary only","Additive only","Neither/NS")))
pal <- c("Both sig"="#2ca02c","Primary only"="#d62728","Additive only"="#ff7f0e","Neither/NS"="grey80")
p_sc <- ggplot(m, aes(x = logFC_primary, y = logFC_add, colour = cat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey85") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey85") +
  geom_point(data = subset(m, cat == "Neither/NS"), size = 0.4, alpha = 0.3) +
  geom_point(data = subset(m, cat != "Neither/NS"), size = 1.1, alpha = 0.7) +
  scale_colour_manual(values = pal) +
  labs(title = "BSJ logFC: primary vs additive (FinalAnalysis 01)",
       subtitle = sprintf("rho=%.3f overall, %.3f among primary sig", rho_all, rho_sig),
       x = "primary logFC (~group)", y = "additive logFC (~group + PSD)",
       colour = "") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "scatter_BSJ_primary_vs_additive.pdf"), p_sc, width = 8, height = 6.5)

cat("\nDone — BSJ DE (primary + additive).\n")
