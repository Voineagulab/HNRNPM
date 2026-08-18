suppressPackageStartupMessages({
  library(limma); library(edgeR); library(ggplot2); library(patchwork)
})

# FinalAnalysis item 01 (FSJ): forward-splice junction DE on CIRIquant FSJ counts.
# Universe: starts from the 5,820 BSJ-filtered universe, then additionally
# filters by `>=2 FSJ reads in >=2 samples`. The dropped 270 circRNAs (~4.6%)
# are "extreme linear-collapse / fully circular" cases — no FSJ at the BSJ
# junction in essentially any sample — that cannot be tested for FSJ DE; they
# are flagged separately and inherit a "linear-collapse" interpretation in the
# downstream bucket scheme.
#
# Two models fitted in parallel:
#   primary  : ~ group         (publication-primary for the FSJ axis)
#   additive : ~ group + PSD   (timepoint sensitivity comparator)

ROOT     <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
COUNTS_F <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/fsj_matrix_CIRIquant_HnrnpM_PSD58.tsv")
MULTIQC  <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/NF_LAUNCH_circrna_230420/results/multiqc/multiqc_data/multiqc_general_stats.txt"
BSJ_F    <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv")
OUT_DIR  <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ANNOT_COLS <- c("gene_symbol","chr","start","end","strand","circ_type","gene_id")

raw <- read.table(COUNTS_F, header = TRUE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE)
annot   <- raw[, c("circRNA", ANNOT_COLS), drop = FALSE]
counts  <- round(as.matrix(raw[, !(colnames(raw) %in% c("circRNA", ANNOT_COLS))]))
rownames(counts) <- raw$circRNA
samples <- colnames(counts)
cat("Loaded FSJ counts (BSJ-filtered universe):", nrow(counts),
    "circRNAs x", ncol(counts), "samples\n")

# Apply FSJ abundance filter
keep_mask <- rowSums(counts >= 2) >= 2
n_pre <- nrow(counts); n_post <- sum(keep_mask)
cat(sprintf("FSJ abundance filter (>=2 reads in >=2 samples): %d / %d (%.1f%%)\n",
            n_post, n_pre, 100*n_post/n_pre))

# Identify and save the "fail" set (extreme linear-collapse / fully-circular)
fail_circs <- rownames(counts)[!keep_mask]
fail_tab <- annot[match(fail_circs, annot$circRNA), ]
fail_tab$fsj_total_7samples <- rowSums(counts[fail_circs, , drop = FALSE])
fail_tab$fsj_max_per_sample <- apply(counts[fail_circs, , drop = FALSE], 1, max)
fail_tab$n_zero_samples     <- rowSums(counts[fail_circs, , drop = FALSE] == 0)
write.table(fail_tab[order(fail_tab$fsj_total_7samples), ],
            file.path(OUT_DIR, "fsj_filter_fail_circRNAs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  -> wrote %d circRNAs to fsj_filter_fail_circRNAs.tsv (linear-collapse candidates)\n",
            length(fail_circs)))

counts <- counts[keep_mask, , drop = FALSE]
annot  <- annot[match(rownames(counts), annot$circRNA), , drop = FALSE]

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
dge <- DGEList(counts = counts, lib.size = lib_sizes, group = group)

# CPM table
cpm_mat <- cpm(dge)
cpm_out <- cbind(circRNA = rownames(cpm_mat),
                 annot[match(rownames(cpm_mat), annot$circRNA), ANNOT_COLS],
                 as.data.frame(cpm_mat))
write.table(cpm_out, file.path(OUT_DIR, "cpm_FSJ_HnrnpM_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

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
res_pri <- fit_de(design_pri, "groupgHnrnpM", "FSJ PRIMARY ~group")
write.table(res_pri, file.path(OUT_DIR, "limma_FSJ_HnrnpM_vs_NEG_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Additive: ~ group + PSD ----
design_add <- model.matrix(~ group + PSD); rownames(design_add) <- samples
res_add <- fit_de(design_add, "groupgHnrnpM", "FSJ ADDITIVE ~group+PSD")
write.table(res_add, file.path(OUT_DIR, "limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# PSD coefficient table from the additive model
v_add   <- voom(dge, design_add, plot = FALSE)
fit_add <- eBayes(lmFit(v_add, design_add))
res_PSD <- topTable(fit_add, coef = "PSD8", number = Inf, sort.by = "P")
res_PSD_out <- cbind(circRNA = rownames(res_PSD),
                     annot[match(rownames(res_PSD), annot$circRNA), ANNOT_COLS],
                     res_PSD)
write.table(res_PSD_out, file.path(OUT_DIR, "limma_FSJ_PSDcoef_PSD58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
n_psd <- sum(res_PSD$adj.P.Val < 0.05, na.rm = TRUE)
cat(sprintf("[PSD coef] PSD8 vs PSD5: sig at adj.P<0.05 = %d\n", n_psd))

# ---- Primary vs additive overlap ----
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
cat(sprintf("\n--- FSJ primary vs additive overlap ---\nprimary sig: %d\nadditive sig: %d\nboth (same dir): %d\nprimary-only: %d\nadditive-only: %d\nflipped: %d\nrho (all): %.3f\nrho (primary sig): %.3f\n",
            sum(m$sig_pri, na.rm=TRUE), sum(m$sig_add, na.rm=TRUE),
            both, pri_only, add_only, flip, rho_all, rho_sig))
ov <- data.frame(
  metric = c("FSJ_primary_sig","FSJ_additive_sig","both_sig_same_dir","primary_only","additive_only",
             "sig_both_flipped","spearman_logFC_all","spearman_logFC_primary_sig",
             "universe_after_FSJ_filter","FSJ_PSD_coef_sig","FSJ_filter_fail_count"),
  value  = c(sum(m$sig_pri, na.rm=TRUE), sum(m$sig_add, na.rm=TRUE),
             both, pri_only, add_only, flip,
             round(rho_all, 4), round(rho_sig, 4),
             nrow(res_pri), n_psd, length(fail_circs)))
write.table(ov, file.path(OUT_DIR, "FSJ_primary_vs_additive_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Compare to BSJ DE: cross-tabulate ----
bsj <- read.table(BSJ_F, sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
cross <- merge(bsj[, c("circRNA","gene_symbol","logFC","P.Value","adj.P.Val")],
               res_pri[, c("circRNA","logFC","P.Value","adj.P.Val")],
               by = "circRNA", suffixes = c("_BSJ","_FSJ"))
# BSJ-sig: adj.P < 0.05 (strict). FSJ-sig: SOFT rule — raw P < 0.05 AND
# |logFC_FSJ| >= log2(1.5). The FSJ per-axis FDR is conservative for low-
# coverage events; the soft rule captures biologically meaningful FSJ shifts
# while still requiring statistical evidence.
LFC_T_SOFT <- log2(1.5)
cross$sig_BSJ <- cross$adj.P.Val_BSJ < 0.05
cross$sig_FSJ <- (cross$P.Value_FSJ < 0.05) & (abs(cross$logFC_FSJ) >= LFC_T_SOFT)
cross$dir_BSJ <- sign(cross$logFC_BSJ)
cross$dir_FSJ <- sign(cross$logFC_FSJ)

# Mechanistic class — 4-way under BSJ-DE only
cross$class <- with(cross, ifelse(
  sig_BSJ & !sig_FSJ, "independent_backsplicing",
  ifelse(sig_BSJ & sig_FSJ & dir_BSJ == dir_FSJ, "co_regulated",
  ifelse(sig_BSJ & sig_FSJ & dir_BSJ != dir_FSJ, "opposite_regulation",
  ifelse(!sig_BSJ & sig_FSJ, "FSJ_only_NS_BSJ", "neither")))))

# Among BSJ-DE circRNAs only
bsj_sig <- cross[cross$sig_BSJ, ]
cat(sprintf("\n--- BSJ-DE circRNAs cross-classified by FSJ status (n=%d in cross-table) ---\n", nrow(bsj_sig)))
print(table(bsj_sig$class))
write.table(cross, file.path(OUT_DIR, "BSJ_FSJ_classification.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# How many BSJ-DE circRNAs are missing from FSJ test (FSJ-filter-fail)?
bsj_sig_ids <- bsj$circRNA[bsj$adj.P.Val < 0.05]
in_fsj_test <- sum(bsj_sig_ids %in% rownames(counts))
missing <- length(bsj_sig_ids) - in_fsj_test
cat(sprintf("\nBSJ-DE circRNAs: %d total; %d in FSJ-test universe; %d FSJ-filter-fail (auto-class as linear-collapse)\n",
            length(bsj_sig_ids), in_fsj_test, missing))

# ---- Volcano + scatter ----
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
ggsave(file.path(OUT_DIR, "volcano_FSJ_HnrnpM_PSD58.pdf"),         volcano(res_pri, "limma-voom FSJ: ~group"),         width = 8, height = 6)
ggsave(file.path(OUT_DIR, "volcano_FSJ_HnrnpM_PSD58_PSDadd.pdf"), volcano(res_add, "limma-voom FSJ: ~group + PSD"), width = 8, height = 6)

# BSJ vs FSJ logFC scatter, BSJ-DE-coloured
plot_df <- cross[!is.na(cross$logFC_BSJ) & !is.na(cross$logFC_FSJ), ]
plot_df$facet_cat <- factor(ifelse(plot_df$sig_BSJ, plot_df$class, "BSJ_NS"),
                            levels = c("independent_backsplicing","co_regulated","opposite_regulation","FSJ_only_NS_BSJ","neither","BSJ_NS"))
pal2 <- c("independent_backsplicing"="#1f77b4","co_regulated"="#2ca02c",
          "opposite_regulation"="#d62728","FSJ_only_NS_BSJ"="#9467bd",
          "neither"="grey80","BSJ_NS"="grey85")
p_bsj_fsj <- ggplot(plot_df, aes(x = logFC_BSJ, y = logFC_FSJ, colour = facet_cat)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = subset(plot_df, !sig_BSJ), size = 0.4, alpha = 0.25) +
  geom_point(data = subset(plot_df, sig_BSJ), size = 1.2, alpha = 0.8) +
  scale_colour_manual(values = pal2) +
  labs(title = "BSJ vs FSJ logFC (~group)",
       subtitle = "BSJ-DE circRNAs coloured by mechanistic class",
       x = "BSJ logFC", y = "FSJ logFC", colour = "") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "scatter_BSJ_vs_FSJ_logFC.pdf"), p_bsj_fsj, width = 9, height = 6.5)

cat("\nDone — FSJ DE.\n")
