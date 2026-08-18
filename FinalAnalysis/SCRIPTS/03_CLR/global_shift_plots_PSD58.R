suppressPackageStartupMessages({
  library(ggplot2); library(tidyr); library(dplyr); library(patchwork)
  library(RColorBrewer); library(ggrepel); library(tibble)
})

# FinalAnalysis item 03: visualise the GLOBAL SHIFT in circRNA formation
# (circular-to-linear ratio, CLR) between NEG and HnrnpM.
#
# Three plots in one combined PDF:
#   A. Per-circRNA mean CLR scatter (NEG vs HnrnpM, y=x diagonal annotated)
#   B. CLR density by group + by sample (detected only, log-x)
#   C. CLR distribution by 4 bins (<0.1 / 0.1-0.5 / 0.5-0.9 / >0.9)

ROOT    <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS"
CLR_F   <- file.path(ROOT, "03_CLR", "clr_matrix_HnrnpM_PSD58.tsv")
BSJ_DE  <- file.path(ROOT, "01_circRNA_DE", "limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv")
OUT_DIR <- file.path(ROOT, "03_CLR")

dat <- read.table(CLR_F, sep = "\t", header = TRUE, check.names = FALSE)
meta_cols <- c("circRNA","gene_symbol","chr","start","end","strand","circ_type","gene_id")
mat <- as.matrix(dat[, setdiff(colnames(dat), meta_cols)])
rownames(mat) <- dat$circRNA

HM  <- c("5_gHnrnpM_PSD5_S20","6_gHnrnpM_PSD5_S21","7_gHnrnpM_PSD8_S22")
NEG <- c("9_gNEG4_PSD5_S24","11_gNEG4_PSD8_S3","12_gNEG4_PSD5_S4","13_gNEG4_PSD8_S5")
sample_order <- c(NEG, HM)

# BSJ-DE circRNAs from step 01 (used to colour the scatter)
bsj_de <- read.table(BSJ_DE, sep = "\t", header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE)
sig_bsj <- bsj_de$circRNA[bsj_de$adj.P.Val < 0.05]

# ---- (A) Per-circRNA mean CLR scatter ----
scat <- data.frame(
  circRNA    = rownames(mat),
  gene_symbol= dat$gene_symbol,
  mean_NEG   = rowMeans(mat[, NEG]),
  mean_HM    = rowMeans(mat[, HM]),
  stringsAsFactors = FALSE
)
scat$sig <- scat$circRNA %in% sig_bsj
scat <- scat[order(scat$sig), ]
n_above <- sum(scat$mean_HM > scat$mean_NEG)
pct_above <- round(100 * n_above / nrow(scat), 1)

p_scatter <- ggplot(scat, aes(x = mean_NEG, y = mean_HM)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = sig, size = sig, alpha = sig)) +
  geom_text_repel(data = subset(scat, sig),
                  aes(label = gene_symbol),
                  size = 2.7, min.segment.length = 0, max.overlaps = 25, seed = 1) +
  scale_colour_manual(values = c("FALSE" = "grey75", "TRUE" = "#d73027"),
                      labels = c("FALSE"="NS", "TRUE"="BSJ-DE (adj.P<0.05)"),
                      name = "") +
  scale_size_manual(values = c("FALSE"=0.6, "TRUE"=2.2), guide = "none") +
  scale_alpha_manual(values = c("FALSE"=0.35, "TRUE"=0.95), guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "A. Mean CLR per circRNA: NEG4 vs HnrnpM",
       subtitle = sprintf("%d circRNAs; %d (%.1f%%) above the y=x diagonal",
                          nrow(scat), n_above, pct_above),
       x = "Mean CLR in NEG4 (n=4)",
       y = "Mean CLR in HnrnpM (n=3)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "scatter_mean_CLR_HnrnpM_PSD58.pdf"),
       p_scatter, width = 7, height = 7.4)

# ---- (B) CLR density by group + by sample ----
long <- as.data.frame(mat) %>% tibble::rownames_to_column("circRNA") %>%
  pivot_longer(-circRNA, names_to = "sample", values_to = "CLR") %>%
  mutate(group  = factor(ifelse(grepl("gHnrnpM", sample), "HnrnpM", "NEG4"),
                         levels = c("NEG4","HnrnpM")),
         sample = factor(sample, levels = sample_order)) %>%
  filter(CLR > 0)

group_pal  <- c(NEG4 = "#4C78A8", HnrnpM = "#E45756")
sample_pal <- c(setNames(brewer.pal(9, "Blues")[5:9][seq_along(NEG)], NEG),
                setNames(brewer.pal(7, "Reds")[4:7][seq_along(HM)],   HM))

p_dens_group <- ggplot(long, aes(x = CLR, colour = group, fill = group)) +
  geom_density(alpha = 0.30, linewidth = 0.8) +
  scale_colour_manual(values = group_pal) +
  scale_fill_manual(values = group_pal) +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1), limits = c(0.001, 1)) +
  labs(title = "B. CLR density by group (detected only, log-x)",
       subtitle = sprintf("%d detected (CLR>0) across 5,820 circRNAs x 7 samples", nrow(long)),
       x = "CLR (log scale)", y = "Density", colour = "Group", fill = "Group") +
  theme_bw(base_size = 11) + theme(legend.position = "right")

p_dens_sample <- ggplot(long, aes(x = CLR, colour = sample, linetype = group, group = sample)) +
  geom_density(linewidth = 0.7) +
  scale_colour_manual(values = sample_pal) +
  scale_linetype_manual(values = c(NEG4 = "dashed", HnrnpM = "solid")) +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1), limits = c(0.001, 1)) +
  labs(title = "C. CLR density by sample (log-x)",
       subtitle = "Blues = NEG4 (dashed); Reds = HnrnpM (solid)",
       x = "CLR (log scale)", y = "Density", colour = "Sample", linetype = "Group") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.text = element_text(size = 8))

# ---- (D) 4-bin histogram ----
bin_breaks <- c(0, 0.1, 0.5, 0.9, 1.001)
bin_labels <- c("<0.1", "0.1-0.5", "0.5-0.9", ">0.9")
bin_df <- long %>%
  mutate(bin = cut(CLR, breaks = bin_breaks, labels = bin_labels,
                   include.lowest = TRUE, right = FALSE)) %>%
  count(group, bin, name = "n") %>%
  group_by(group) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()

p_hist <- ggplot(bin_df, aes(x = bin, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "white") +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = group_pal) +
  labs(title = "D. CLR distribution by bin (% of detected, per group)",
       x = "CLR bin", y = "% detected per group", fill = "Group") +
  theme_bw(base_size = 11)

# Save density + histogram in one PDF
combined_dh <- (p_dens_group / p_dens_sample / p_hist)
ggsave(file.path(OUT_DIR, "density_CLR_HnrnpM_PSD58.pdf"),
       combined_dh, width = 9, height = 13)

# Also: a single 4-panel summary PDF combining scatter + density + histogram
combined_all <- (p_scatter | (p_dens_group / p_dens_sample / p_hist)) +
                plot_layout(widths = c(1, 1.1))
ggsave(file.path(OUT_DIR, "global_shift_CLR_combined_PSD58.pdf"),
       combined_all, width = 16, height = 11)

cat("Wrote:\n",
    " -", file.path(OUT_DIR, "scatter_mean_CLR_HnrnpM_PSD58.pdf"), "\n",
    " -", file.path(OUT_DIR, "density_CLR_HnrnpM_PSD58.pdf"), "\n",
    " -", file.path(OUT_DIR, "global_shift_CLR_combined_PSD58.pdf"), "\n")

cat(sprintf("\nDirection summary: %d / %d circRNAs (%.1f%%) above y=x diagonal\n",
            n_above, nrow(scat), pct_above))
cat("\nPer-bin distribution (% of detected, per group):\n")
print(bin_df %>% select(group, bin, n, pct) %>% as.data.frame(), row.names = FALSE)
