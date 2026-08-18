#!/usr/bin/env Rscript
# Compare two FSJ-significance definitions for the BSJ x FSJ mechanistic
# classification of BSJ-DE circRNAs (additive `~ group + PSD` model):
#   * STRICT: FSJ adj.P < 0.05
#   * SOFT  : FSJ raw P < 0.05 AND |FSJ logFC| >= log2(1.5)
#
# Outputs (RESULTS/01_circRNA_DE/FSJ/ThresholdSelection/):
#   FSJ_threshold_summary.tsv          — # sig FSJ events + class counts
#   BSJ_DE_dual_classification.tsv     — per BSJ-DE circRNA, strict + soft class
#   boxplots_BSJ_DE_threshold_comparison.pdf  — 56-panel BSJ/FSJ CPM boxplots

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

ROOT <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE"
FSJ  <- file.path(ROOT, "FSJ")
OUT  <- file.path(FSJ,  "ThresholdSelection")

bsj_de   <- read.table(file.path(ROOT, "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"), sep="\t", header=TRUE, check.names=FALSE)
fsj_de   <- read.table(file.path(FSJ,  "limma_FSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv"), sep="\t", header=TRUE, check.names=FALSE)
fsj_fail <- read.table(file.path(FSJ,  "fsj_filter_fail_circRNAs.tsv"),             sep="\t", header=TRUE, check.names=FALSE)
fsj_fail_set <- fsj_fail$circRNA

bsj_cpm  <- read.table(file.path(ROOT, "cpm_BSJ_HnrnpM_PSD58.tsv"),                 sep="\t", header=TRUE, check.names=FALSE)
fsj_cpm  <- read.table(file.path(FSJ,  "cpm_FSJ_HnrnpM_PSD58.tsv"),                 sep="\t", header=TRUE, check.names=FALSE)

LFC_T_SOFT <- log2(1.5)

# ----- universe-level FSJ significance counts ---------------------------------
n_fsj_tested <- nrow(fsj_de)
n_fsj_strict <- sum(fsj_de$adj.P.Val < 0.05, na.rm=TRUE)
n_fsj_soft   <- sum(fsj_de$P.Value   < 0.05 & abs(fsj_de$logFC) >= LFC_T_SOFT, na.rm=TRUE)

# ----- merge BSJ-DE with FSJ stats -------------------------------------------
bsj_keep <- bsj_de[, c("circRNA","gene_symbol","logFC","adj.P.Val")]
names(bsj_keep) <- c("circRNA","gene_symbol","logFC_BSJ","adjP_BSJ")
fsj_keep <- fsj_de[, c("circRNA","logFC","P.Value","adj.P.Val")]
names(fsj_keep) <- c("circRNA","logFC_FSJ","rawP_FSJ","adjP_FSJ")

m <- merge(bsj_keep, fsj_keep, by="circRNA", all.x=TRUE)
m$sig_BSJ          <- m$adjP_BSJ < 0.05
m$in_fsj_test      <- !is.na(m$adjP_FSJ)
m$fsj_filter_fail  <- m$circRNA %in% fsj_fail_set

m$sig_FSJ_strict <- !is.na(m$adjP_FSJ) & m$adjP_FSJ < 0.05
m$sig_FSJ_soft   <- !is.na(m$rawP_FSJ) & m$rawP_FSJ < 0.05 &
                     !is.na(m$logFC_FSJ) & abs(m$logFC_FSJ) >= LFC_T_SOFT

classify <- function(sig_BSJ, sig_FSJ, dir_BSJ, dir_FSJ, in_fsj, fail) {
  if (sig_BSJ && fail)        return("extreme_linear_collapse")
  if (sig_BSJ && !in_fsj)     return("extreme_linear_collapse")
  if (sig_BSJ &&  sig_FSJ) {
    if (dir_BSJ == dir_FSJ)   return("co_regulated")
    else                      return("opposite_regulation")
  }
  if (sig_BSJ && !sig_FSJ)    return("independent_backsplicing")
  if (!sig_BSJ && sig_FSJ)    return("FSJ_only_NS_BSJ")
  return("neither")
}

m$dir_BSJ <- sign(m$logFC_BSJ)
m$dir_FSJ <- sign(m$logFC_FSJ)

m$class_strict <- mapply(classify, m$sig_BSJ, m$sig_FSJ_strict, m$dir_BSJ, m$dir_FSJ,
                          m$in_fsj_test, m$fsj_filter_fail)
m$class_soft   <- mapply(classify, m$sig_BSJ, m$sig_FSJ_soft,   m$dir_BSJ, m$dir_FSJ,
                          m$in_fsj_test, m$fsj_filter_fail)

# ----- subset to BSJ-DE -------------------------------------------------------
bsj_sig <- m[m$sig_BSJ, ]
class_levels <- c("independent_backsplicing","co_regulated","opposite_regulation","extreme_linear_collapse")
strict_counts <- table(factor(bsj_sig$class_strict, levels=class_levels))
soft_counts   <- table(factor(bsj_sig$class_soft,   levels=class_levels))

# ----- summary TSV -----------------------------------------------------------
summary_tbl <- data.frame(
  metric = c("FSJ events tested",
             "FSJ events significant",
             paste0("BSJ-DE: ", class_levels)),
  strict_adjP_lt_0.05 = c(n_fsj_tested, n_fsj_strict,
                          as.integer(strict_counts)),
  soft_rawP_and_FC = c(n_fsj_tested, n_fsj_soft,
                       as.integer(soft_counts))
)
write.table(summary_tbl, file.path(OUT, "FSJ_threshold_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote:", file.path(OUT, "FSJ_threshold_summary.tsv"), "\n")
print(summary_tbl, row.names=FALSE)

# ----- per-circRNA dual classification ---------------------------------------
out_dual <- bsj_sig[, c("circRNA","gene_symbol","logFC_BSJ","adjP_BSJ",
                         "logFC_FSJ","rawP_FSJ","adjP_FSJ",
                         "in_fsj_test","fsj_filter_fail",
                         "sig_FSJ_strict","sig_FSJ_soft",
                         "class_strict","class_soft")]
out_dual$class_changes <- out_dual$class_strict != out_dual$class_soft
out_dual <- out_dual[order(out_dual$class_soft, out_dual$adjP_BSJ), ]
write.table(out_dual, file.path(OUT, "BSJ_DE_dual_classification.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWrote:", file.path(OUT, "BSJ_DE_dual_classification.tsv"), "\n")
cat("Class changes between strict and soft:", sum(out_dual$class_changes), "/", nrow(out_dual), "\n")

# ----- box-plots -------------------------------------------------------------
neg_cols <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
hm_cols  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")

# Build per-circRNA tidy data for the 56 BSJ-DE circRNAs.
mk_long <- function(cpm_df, junction, circs) {
  d <- cpm_df[cpm_df$circRNA %in% circs, c("circRNA","gene_symbol", neg_cols, hm_cols)]
  long <- pivot_longer(d, all_of(c(neg_cols, hm_cols)),
                       names_to="sample", values_to="cpm")
  long$group <- ifelse(long$sample %in% hm_cols, "KO", "NEG")
  long$junction <- junction
  long$log2cpm  <- log2(long$cpm + 1)
  long
}

circs <- bsj_sig$circRNA
bsj_long <- mk_long(bsj_cpm, "BSJ", circs)
fsj_long <- mk_long(fsj_cpm, "FSJ", circs)
plot_df  <- bind_rows(bsj_long, fsj_long)
plot_df$jg <- factor(paste(plot_df$junction, plot_df$group, sep="-"),
                     levels=c("BSJ-NEG","BSJ-KO","FSJ-NEG","FSJ-KO"))

# Per-panel title encoding both classifications.
abbr <- c(independent_backsplicing="indep",
          co_regulated="co-reg",
          opposite_regulation="opp",
          extreme_linear_collapse="collapse",
          FSJ_only_NS_BSJ="fsj-only",
          neither="neither")
ord <- bsj_sig %>%
  arrange(class_soft, adjP_BSJ) %>%
  group_by(gene_symbol) %>%
  mutate(gene_label = if (n() > 1) paste0(gene_symbol, "#", row_number()) else gene_symbol) %>%
  ungroup() %>%
  mutate(panel_title = sprintf("%s\nstrict: %s | soft: %s",
                                gene_label,
                                abbr[class_strict], abbr[class_soft]))
plot_df <- plot_df %>%
  left_join(ord[, c("circRNA","panel_title","class_soft")], by="circRNA")
plot_df$panel_title <- factor(plot_df$panel_title, levels=ord$panel_title)

pal <- c("BSJ-NEG"="#4C78A8","BSJ-KO"="#E45756",
         "FSJ-NEG"="#9ec3e3","FSJ-KO"="#f5a8a8")

p <- ggplot(plot_df, aes(x=jg, y=log2cpm, fill=jg)) +
  geom_boxplot(outlier.shape=NA, linewidth=0.3, colour="grey20") +
  geom_jitter(width=0.15, size=0.7, alpha=0.8, colour="grey20") +
  scale_fill_manual(values=pal, name=NULL) +
  facet_wrap(~ panel_title, ncol=7, scales="free_y") +
  labs(title="BSJ-DE circRNAs (additive ~group+PSD): BSJ vs FSJ CPM",
       subtitle="strict = FSJ adj.P<0.05 | soft = FSJ raw P<0.05 AND |FSJ logFC|>=log2(1.5)",
       x=NULL, y="log2(CPM + 1)") +
  theme_bw(base_size=8) +
  theme(strip.text=element_text(size=6.2, lineheight=0.9),
        axis.text.x=element_text(angle=45, hjust=1, size=6),
        legend.position="bottom",
        panel.spacing=unit(0.4, "lines"))

ggsave(file.path(OUT, "boxplots_BSJ_DE_threshold_comparison.pdf"),
       p, width=14, height=18)
cat("\nWrote:", file.path(OUT, "boxplots_BSJ_DE_threshold_comparison.pdf"), "\n")
