#!/usr/bin/env Rscript
# =============================================================================
# d2_CLR_violins.R  (log10 version)
# -----------------------------------------------------------------------------
# Two CLR violin plots for the manuscript, built from the FinalAnalysis CLR
# matrix (RESULTS/03_CLR/clr_matrix_HnrnpM_PSD58.tsv). Read-only on all inputs;
# all outputs are written into this folder (d2_CLR) only.
#
#   Plot 1  CLR per sample, 7 violins (4 control, 3 KD), with a horizontal line
#           at each sample's mean CLR.
#   Plot 2  Per-circRNA mean CLR in the control and KD groups, 2 violins (PAIRED).
#
# LOG TRANSFORM NOTE
#   The y axis is log10(CLR). Because CLR = 0 (undetected) cannot be shown on a
#   log axis, the VIOLINS show detected values only (CLR > 0). The MEAN LINE is
#   deliberately kept as the arithmetic mean over ALL 5,820 circRNAs (including
#   zeros), because that is the exact quantity used by the per-sample Welch t /
#   Mann-Whitney tests (Plot 1) and the per-circRNA paired Wilcoxon / sign tests
#   (Plot 2) in global_CLR_shift_tests.tsv. So the violin shows the shape of the
#   detected distribution on a log scale, while the line marks the tested mean.
#
# Stats annotated are read from global_CLR_shift_tests.tsv so they stay faithful.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
})

# ---- paths (Docker / RStudio mount) ----
ROOT  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS"
CLR_F <- file.path(ROOT, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
TESTS <- file.path(ROOT, "03_CLR", "global_CLR_shift_tests.tsv")
OUT   <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d2_CLR"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- sample groups (control first, then KD) ----
HM  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
NEG <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
sample_order <- c(NEG, HM)
meta_cols <- c("circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id")

relabel <- function(x) gsub("gHnrnpM","gHNRNPM", gsub("gNEG4","gCTRL", x))
group_pal <- c(gCTRL = "#4C78A8", gHNRNPM = "#E45756")
log_breaks <- c(0.001, 0.01, 0.1, 1)

# ---- load CLR matrix ----
dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE)
samples <- setdiff(colnames(dat), meta_cols)
mat <- as.matrix(dat[, samples]); rownames(mat) <- dat$circRNA
cat(sprintf("CLR matrix: %d circRNAs x %d samples\n", nrow(mat), ncol(mat)))

# ---- read the tests table for faithful annotation ----
tests <- read.table(TESTS, sep = "\t", header = TRUE, check.names = FALSE)
getp <- function(pat) {
  hit <- tests$p[grepl(pat, tests$test)]
  if (length(hit) == 0) NA_real_ else as.numeric(hit[1])
}
p_welch  <- getp("Welch_t")
p_mwu    <- getp("Mann-Whitney")
p_wilcox <- getp("Wilcoxon")
p_sign   <- getp("sign_test")

fp <- function(p) {                                   # pretty p-value
  if (is.na(p)) return("NA")
  if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "g", digits = 2)
}

# =============================================================================
# PLOT 1 — CLR per sample (7 violins, log10 y) + per-sample mean line
# =============================================================================
long1 <- as.data.frame(mat) |>
  tibble::rownames_to_column("circRNA") |>
  pivot_longer(-circRNA, names_to = "sample", values_to = "CLR") |>
  mutate(group = ifelse(grepl("gHnrnpM", sample), "gHNRNPM", "gCTRL"),
         sample_lab = factor(relabel(sample), levels = relabel(sample_order)))

# violins: detected only (CLR > 0), because zeros cannot be shown on a log axis
long1_det <- dplyr::filter(long1, CLR > 0)
cat(sprintf("Plot 1: %d of %d (circRNA x sample) entries are detected (CLR>0) and shown in the violins; %d zeros omitted by the log axis.\n",
            nrow(long1_det), nrow(long1), nrow(long1) - nrow(long1_det)))

# mean line: arithmetic mean over ALL circRNAs (incl zeros) = the tested quantity
samp_mean <- long1 |> group_by(sample_lab, group) |>
  summarise(mean_CLR = mean(CLR), .groups = "drop") |>
  mutate(x = as.integer(sample_lab))
cat("\nPer-sample mean CLR (over all 5,820 circRNAs, incl. zeros):\n")
print(as.data.frame(samp_mean[, c("sample_lab","group","mean_CLR")]), row.names = FALSE)

p1 <- ggplot(long1_det, aes(x = sample_lab, y = CLR, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey30", linewidth = 0.3, alpha = 0.6) +
  geom_segment(data = samp_mean,
               aes(x = x - 0.45, xend = x + 0.45, y = mean_CLR, yend = mean_CLR),
               inherit.aes = FALSE, linewidth = 0.8, colour = "black") +
  geom_text(data = samp_mean,
            aes(x = x, y = mean_CLR, label = sprintf("%.3f", mean_CLR)),
            inherit.aes = FALSE, vjust = -0.6, size = 2.6) +
  scale_fill_manual(values = group_pal, name = "Group") +
  scale_y_log10(breaks = log_breaks, labels = log_breaks) +
  labs(title = "Circular-to-linear ratio (CLR) per sample (log10 scale)",
       subtitle = sprintf("Violins = detected CLR (>0). Black line = per-sample mean CLR over all 5,820 circRNAs. Group difference (unit = sample): Welch t P = %s; Mann-Whitney U P = %s (3 KD vs 4 control).",
                          fp(p_welch), fp(p_mwu)),
       x = NULL, y = "CLR (log10 scale)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.subtitle = element_text(size = 6.5))

ggsave(file.path(OUT, "Fig_CLR_per_sample_violin_log10.pdf"), p1, width = 18, height = 9, units = "cm", device = cairo_pdf)
ggsave(file.path(OUT, "Fig_CLR_per_sample_violin_log10.png"), p1, width = 18, height = 9, units = "cm", dpi = 200)

# =============================================================================
# PLOT 2 — per-circRNA mean CLR by group (2 violins, log10 y, paired)
# =============================================================================
mean_NEG <- rowMeans(mat[, NEG]); mean_HM <- rowMeans(mat[, HM])
long2 <- tibble(circRNA = rownames(mat),
                gCTRL = mean_NEG, gHNRNPM = mean_HM) |>
  pivot_longer(c(gCTRL, gHNRNPM), names_to = "group", values_to = "mean_CLR") |>
  mutate(group = factor(group, levels = c("gCTRL","gHNRNPM")))

# violins: detected only (per-circRNA group mean > 0)
long2_det <- dplyr::filter(long2, mean_CLR > 0)
cat(sprintf("\nPlot 2: %d of %d per-circRNA group means are > 0 and shown; %d zeros omitted by the log axis.\n",
            nrow(long2_det), nrow(long2), nrow(long2) - nrow(long2_det)))

# mean line: arithmetic group mean over ALL circRNAs (incl zeros) = tested quantity
grp_stat <- long2 |> group_by(group) |>
  summarise(mean = mean(mean_CLR), .groups = "drop") |>
  mutate(x = as.integer(group))
pct_above <- 100 * mean(mean_HM > mean_NEG)
cat(sprintf("Group-level mean of per-circRNA mean CLR: gCTRL %.4f, gHNRNPM %.4f\n",
            grp_stat$mean[grp_stat$group=="gCTRL"], grp_stat$mean[grp_stat$group=="gHNRNPM"]))
cat(sprintf("Per-circRNA: mean_HM > mean_NEG in %.1f%% of circRNAs\n", pct_above))

p2 <- ggplot(long2_det, aes(x = group, y = mean_CLR, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey30", linewidth = 0.3, alpha = 0.6) +
  geom_segment(data = grp_stat,
               aes(x = x - 0.4, xend = x + 0.4, y = mean, yend = mean),
               inherit.aes = FALSE, linewidth = 0.8, colour = "black") +
  geom_text(data = grp_stat,
            aes(x = x, y = mean, label = sprintf("mean %.3f", mean)),
            inherit.aes = FALSE, vjust = -0.6, size = 2.6) +
  scale_fill_manual(values = group_pal, name = "Group", guide = "none") +
  scale_y_log10(breaks = log_breaks, labels = log_breaks) +
  labs(title = "Per-circRNA mean CLR by group (log10 scale)",
       subtitle = sprintf("Violins = per-circRNA mean CLR (>0); black line = arithmetic group mean over all circRNAs. PAIRED comparison: Wilcoxon signed-rank P = %s; %.1f%% of circRNAs higher in KD (sign test P = %s).",
                          fp(p_wilcox), pct_above, fp(p_sign)),
       x = NULL, y = "Per-circRNA mean CLR (log10 scale)") +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 6.5))

ggsave(file.path(OUT, "Fig_CLR_per_group_violin_log10.pdf"), p2, width = 10, height = 9, units = "cm", device = cairo_pdf)
ggsave(file.path(OUT, "Fig_CLR_per_group_violin_log10.png"), p2, width = 10, height = 9, units = "cm", dpi = 200)

cat("\nWrote to:", OUT, "\n")
cat(list.files(OUT, pattern = "log10\\.(pdf|png)$"), sep = "\n")
