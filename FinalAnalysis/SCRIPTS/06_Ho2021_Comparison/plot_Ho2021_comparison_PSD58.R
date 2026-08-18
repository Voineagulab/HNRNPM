suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})

ROOT <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/06_Ho2021_Comparison"
PG   <- file.path(ROOT, "Ho2021_vs_ours_HnrnpM_DE_PSD58.tsv")
SUM  <- file.path(ROOT, "Ho2021_vs_ours_HnrnpM_summary_PSD58.tsv")
BV   <- file.path(ROOT, "bound_vs_unbound_replication_PSD58.tsv")

m  <- read.table(PG,  sep = "\t", header = TRUE, check.names = FALSE)
s  <- read.table(SUM, sep = "\t", header = TRUE, check.names = FALSE)
bv <- read.table(BV,  sep = "\t", header = TRUE, check.names = FALSE)

# Scatter: Ho mean log2FC vs our limma logFC
m$category <- with(m, ifelse(
  is.na(limma_adj_P), "Not in limma test",
  ifelse(Ho_sig_strict == "True" | Ho_sig_strict == TRUE,
    ifelse(limma_sig == "True" | limma_sig == TRUE,
      ifelse(sign(Ho_mean_log2FC) == sign(limma_logFC), "Concordant (both sig)", "Discordant (both sig)"),
      "Ho-only sig"),
    ifelse(limma_sig == "True" | limma_sig == TRUE, "Our-only sig", "Neither sig"))))
m$category <- factor(m$category, levels = c("Concordant (both sig)","Discordant (both sig)",
                                              "Ho-only sig","Our-only sig","Neither sig","Not in limma test"))
pal <- c("Concordant (both sig)"="#2ca02c","Discordant (both sig)"="#d62728",
         "Ho-only sig"="#1f77b4","Our-only sig"="#ff7f0e",
         "Neither sig"="grey80","Not in limma test"="grey92")
plot_df <- m[!is.na(m$Ho_mean_log2FC) & !is.na(m$limma_logFC), ]
plot_df <- plot_df[order(plot_df$category, decreasing = TRUE), ]

p_scatter <- ggplot(plot_df, aes(x = Ho_mean_log2FC, y = limma_logFC, colour = category)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(data = subset(plot_df, category %in% c("Neither sig","Not in limma test")),
             size = 0.3, alpha = 0.25) +
  geom_point(data = subset(plot_df, !category %in% c("Neither sig","Not in limma test")),
             size = 0.8, alpha = 0.7) +
  scale_colour_manual(values = pal) +
  labs(title = "Ho 2021 vs FinalAnalysis additive (~group + PSD): per-gene logFC",
       subtitle = sprintf("N=%d  ·  Concordant overlap k=%d UP / %d DOWN  ·  P_combined = %.2e (Fisher's method)",
                          s$N, s$k_up, s$k_down, s$P_combined),
       x = "Ho 2021 mean log2FC (B7 + B9)", y = "Our limma logFC (additive)",
       colour = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# Bar plot: bound vs unbound percent concordant
bv$pct_concordant <- as.numeric(bv$pct_concordant)
bv$pct_discordant <- 100 * bv$n_discordant / bv$n_total
bv$stratum <- factor(bv$stratum, levels = c("Bound","Unbound"))
bar_df <- rbind(
  data.frame(stratum = bv$stratum, kind = "Concordant", pct = bv$pct_concordant),
  data.frame(stratum = bv$stratum, kind = "Discordant", pct = bv$pct_discordant)
)
p_bar <- ggplot(bar_df, aes(x = stratum, y = pct, fill = kind)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "white") +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Concordant" = "#2ca02c", "Discordant" = "#d62728")) +
  labs(title = "Ho-strict-sig genes: replication % in our additive DE (Bound vs Unbound)",
       x = NULL, y = "% of Ho strict-sig genes (per stratum)", fill = "") +
  theme_bw(base_size = 11)

combined <- (p_scatter / p_bar) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(title = "FinalAnalysis 06: Ho et al. 2021 cross-cell-type comparison",
                  subtitle = "Anchored on our additive (~ group + PSD) gene-level DE.",
                  theme = theme(plot.title = element_text(size = 12, face = "bold"),
                                plot.subtitle = element_text(size = 10, colour = "grey40")))
ggsave(file.path(ROOT, "Ho2021_vs_ours_HnrnpM_combined_PSD58.pdf"), combined, width = 10, height = 11)
cat("Wrote:", file.path(ROOT, "Ho2021_vs_ours_HnrnpM_combined_PSD58.pdf"), "\n")
