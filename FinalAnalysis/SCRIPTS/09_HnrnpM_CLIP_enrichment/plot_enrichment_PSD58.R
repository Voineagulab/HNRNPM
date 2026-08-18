suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

ROOT  <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS"
SIG   <- file.path(ROOT, "09_HnrnpM_CLIP_enrichment", "clip_signal_per_region.tsv")
SCORE <- file.path(ROOT, "09_HnrnpM_CLIP_enrichment", "circRNA_flank_score_PSD58.tsv")
TESTS <- file.path(ROOT, "09_HnrnpM_CLIP_enrichment", "enrichment_test_summary.tsv")
BSJ_ADD <- file.path(ROOT, "01_circRNA_DE", "limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
OUT_DIR <- file.path(ROOT, "09_HnrnpM_CLIP_enrichment")

sig   <- read.table(SIG,   sep = "\t", header = TRUE, check.names = FALSE)
sc    <- read.table(SCORE, sep = "\t", header = TRUE, check.names = FALSE)
tests <- read.table(TESTS, sep = "\t", header = TRUE, check.names = FALSE)
bsj   <- read.table(BSJ_ADD, sep = "\t", header = TRUE, check.names = FALSE)

# Identify additive BSJ-DE circRNAs (publication primary)
de_circs <- bsj$circRNA[bsj$adj.P.Val < 0.05]
de_n <- sum(sc$circRNA %in% de_circs)
nde_n <- sum(!(sc$circRNA %in% de_circs))

# Panel A: per-region log2(IP/Input) box plot
sig$de_lab <- ifelse(sig$circRNA %in% de_circs, "DE", "non-DE")
sig$panel_label <- with(sig, ifelse(
  region_type == "flank5", paste0(de_lab, " 5' flank"),
  ifelse(region_type == "flank3", paste0(de_lab, " 3' flank"),
  ifelse(region_type == "within_gene_bg", "Within-gene bg",
  ifelse(region_type == "non_host_bg",   "Non-host-gene bg", region_type)))))
sig$panel_label <- factor(sig$panel_label,
  levels = c("DE 5' flank","DE 3' flank","non-DE 5' flank","non-DE 3' flank",
             "Within-gene bg","Non-host-gene bg"))
palette <- c("DE 5' flank" = "#d73027", "DE 3' flank" = "#d73027",
             "non-DE 5' flank" = "#4C78A8", "non-DE 3' flank" = "#4C78A8",
             "Within-gene bg" = "grey50", "Non-host-gene bg" = "grey75")

panel_A <- ggplot(sig, aes(x = panel_label, y = log2_IP_over_Input, fill = panel_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, colour = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = palette, guide = "none") +
  labs(title = "A. log2(IP/Input) per region",
       x = NULL, y = "log2(IP / Input)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel B: per-circRNA mean flank score
sc$is_DE_lab <- ifelse(sc$circRNA %in% de_circs,
                       sprintf("DE (n=%d)", de_n),
                       sprintf("non-DE (n=%d)", nde_n))
p_de_lab <- sprintf("DE (n=%d)", de_n)
p_nde_lab <- sprintf("non-DE (n=%d)", nde_n)
p_anno <- tests$p_one_sided[tests$DE_axis == "BSJ_additive" &
                              tests$test == "1. DE flank vs non-DE flank"]

panel_B <- ggplot(sc, aes(x = is_DE_lab, y = mean_flank_log2, fill = is_DE_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_violin(alpha = 0.5, colour = "grey20", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  scale_fill_manual(values = setNames(c("#d73027","#4C78A8"), c(p_de_lab, p_nde_lab)),
                    guide = "none") +
  labs(title = "B. Per-circRNA mean flank score (DE vs non-DE, additive BSJ-DE)",
       subtitle = sprintf("Wilcoxon rank-sum (one-sided): P = %.2e", p_anno),
       x = NULL, y = "Mean log2(IP/Input) at flank introns") +
  theme_bw(base_size = 11)

combined <- (panel_A / panel_B) +
  plot_annotation(
    title = "HnrnpM eCLIP enrichment at circRNA flanking introns (FinalAnalysis 09)",
    subtitle = "Publication primary: additive BSJ-DE (~ group + PSD).",
    theme = theme(plot.title = element_text(size = 12),
                  plot.subtitle = element_text(size = 9, colour = "grey40"))
  )
ggsave(file.path(OUT_DIR, "enrichment_boxplot_PSD58.pdf"), combined, width = 9, height = 9)
cat("Wrote:", file.path(OUT_DIR, "enrichment_boxplot_PSD58.pdf"), "\n")
