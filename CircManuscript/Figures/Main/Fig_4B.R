# Fig 4B — Ho et al. 2021 cross-cell-type comparison (scatter + bound/unbound bar).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/06_Ho2021_Comparison/plot_Ho2021_comparison_PSD58.R
# Reproduces `Ho2021_vs_ours_HnrnpM_combined_PSD58.pdf`: the per-gene logFC scatter
#   (Ho 2021 vs our additive limma) on top and the bound-vs-unbound replication
#   bar plot below, combined with patchwork. Reads the saved result tables directly.
#
# Strategy (same as Fig_2/Fig_3/Fig_4A): keep ALL original plot parameters, render
#   at the original size, and let ggsave(scale=) squash the figure to 7.75 x 9.562 cm
#   (uniform shrink preserves proportions). Pre-inflate ONLY the axis title/tick
#   fonts + legend so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x
#   larger than nominal; resize it to 7.75 cm wide on the slide (that lands the axis
#   at 8 pt). See Fig_2E_exactsize.R for the born-at-size variant.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES06    <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/06_Ho2021_Comparison")
PG   <- file.path(RES06, "Ho2021_vs_ours_HnrnpM_DE_PSD58.tsv")
SUM  <- file.path(RES06, "Ho2021_vs_ours_HnrnpM_summary_PSD58.tsv")
BV   <- file.path(RES06, "bound_vs_unbound_replication_PSD58.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 9.562     # nominal final PDF size on the slide (70% of the previous 13.66)
ORIG_W_IN <- 10                    # the source rendered the combined figure at width = 10 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Load result tables ----
m  <- read.table(PG,  sep = "\t", header = TRUE, check.names = FALSE)
s  <- read.table(SUM, sep = "\t", header = TRUE, check.names = FALSE)
bv <- read.table(BV,  sep = "\t", header = TRUE, check.names = FALSE)

# ---- Scatter: Ho mean log2FC vs our limma logFC (original parameters preserved) ----
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
       x = "mean log2FC \n(Ho et al, 2021; B7+B9)", y = "limma log2FC \n(This study)",
       colour = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 0.8 * 3)))   # scatter legend dots 3x the plot dots

# ---- Bar plot: bound vs unbound percent concordant ----
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
            position = position_dodge(width = 0.8), vjust = -0.3, size = 4.5) +   # bigger bar data labels
  scale_fill_manual(values = c("Concordant" = "#2ca02c", "Discordant" = "#d62728")) +
  scale_y_continuous(limits = c(0, 35)) +
  labs(title = "Ho-strict-sig genes: replication % in our additive DE (Bound vs Unbound)",
       x = NULL, y = "% of Ho DE genes \n(in both B7 and B9)", fill = "") +
  theme_bw(base_size = 11)

# ---- Combine (scatter over bar); pre-inflate axis/legend on BOTH via `&` ----
combined <- (p_scatter / p_bar) +
  plot_layout(heights = c(1.6, 1.4)) &   # taller bar panel, slightly shorter scatter
  theme(
    axis.title   = element_text(size = 0.75 * axis_render_pt),  # uniform: same size as the bar y-axis label
    axis.text    = element_text(size = 0.75 * axis_render_pt),  # ticks same size as axis labels, both panels
    legend.title = element_text(size = 0.5 * axis_render_pt),   # legend text 50% (-> ~4 pt final)
    legend.text  = element_text(size = 0.5 * axis_render_pt),   # legend text 50% (-> ~4 pt final)
    legend.box.spacing = grid::unit(1, "pt"),                  # tighten gap between plot and legend
    legend.margin      = margin(0, 0, 0, 0),
    plot.margin        = margin(t = 6, r = 2, b = 6, l = 2, unit = "pt")  # vertical space between the two panels
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_4B.pdf"), combined,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_4B.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_4B.pdf)
# ============================================================
# Cross-cell-type comparison of the HNRNPM knockdown gene-expression signature
# against Ho et al. 2021 (PC-3M prostate cancer), anchored on our additive
# (~ group + PSD) gene-level differential expression.
#
# Ho et al. 2021 generated two independent CRISPR HNRNPM-knockdown clones, B7 and
# B9, each compared to a scrambled control. Throughout, a Ho "DE gene" is defined
# by a stringent criterion, significant (q < 0.05) in BOTH the B7 and B9 clones AND
# changing in the same direction in both, rather than in a single clone or a pooled
# test. This yields a high-confidence set of ~2,100 Ho DE genes used for all overlap
# and replication statistics.
#
# Top panel, per-gene scatter of the Ho 2021 mean log2 fold change (average of clones
# B7 and B9) against our limma log2 fold change (additive model). Each point is one
# gene; the dotted lines mark zero on each axis and the dashed diagonal marks equal
# fold change. Points are coloured by significance category (concordant or discordant
# when significant in both datasets, Ho-only, our-only, neither, or not in our test
# set). The subtitle reports the concordant-overlap counts (up and down) and
# P_combined. P_combined is a single overall enrichment p-value obtained by testing
# the up-regulated and down-regulated overlaps SEPARATELY (each a hypergeometric
# enrichment test) and combining the two direction-specific p-values with Fisher's
# method (statistic -2 x [ln P_up + ln P_down], chi-square with 4 degrees of freedom);
# testing the directions separately avoids counting genes that are significant in both
# datasets but move in opposite directions as corroboration.
#
# Bottom panel, replication rate of the Ho DE genes in our additive DE, split by
# whether the gene is HNRNPM-bound (Ho CLIP) versus unbound. Bars show the percentage
# of Ho DE genes in each stratum that replicate (concordant or discordant), with the
# denominator taken within each stratum. Together the panels show that the HNRNPM-KD
# signature is conserved across cell types and that HNRNPM-bound targets replicate at
# a modestly higher rate than unbound genes.
