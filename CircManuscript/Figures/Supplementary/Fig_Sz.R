# Fig Sz — HNRNPM-bound vs unbound replication of Ho 2021 DE genes (one-sided Fisher's exact).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/Other/SCRIPTS/06_Ho2021_Comparison/plot_bound_vs_unbound_PSD58.R
# This is the dedicated graphical counterpart of the Methods sentence
#   "Replication of HNRNPM-bound versus unbound genes was compared by one-sided Fisher's exact test."
#   It reads the same additive-model replication table that Fig_4B uses (FinalAnalysis 06),
#   recomputes the Wilson 95% CIs and the one-sided Fisher's exact test per metric, and
#   annotates the odds ratio + P value on each facet.
#
# Strategy (same as the Fig_2/Fig_3/Fig_4/Fig_5 series): keep ALL original plot parameters
#   (bar/errorbar/text sizes, facets, colours, title/subtitle), render at the original size,
#   and let ggsave(scale=) squash the figure to 7.75 x 6.83 cm (uniform shrink preserves
#   proportions). Pre-inflate ONLY the axis title/tick fonts so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger than
#   nominal; resize it to 7.75 cm wide on the slide. See Fig_2E_exactsize.R.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES06   <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/06_Ho2021_Comparison")
BV      <- file.path(RES06, "bound_vs_unbound_replication_PSD58.tsv")   # same table Fig_4B reads
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Supplementary")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 7                     # the source rendered this figure at width = 7 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- load the additive-model replication strata (Bound vs Unbound) ----
strata <- read.table(BV, sep = "\t", header = TRUE, check.names = FALSE)

# ---- Wilson 95% CI (identical to source) ----
wilson_ci <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA))
  z <- qnorm(1 - (1 - conf) / 2); p <- k / n; denom <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt((p * (1-p) / n) + (z^2 / (4 * n^2))) / denom
  c(c - half, c + half)
}
ci  <- t(mapply(wilson_ci, strata$n_concordant, strata$n_total))
strata$conc_lo <- 100 * ci[,1]; strata$conc_hi <- 100 * ci[,2]
ci2 <- t(mapply(wilson_ci, strata$n_our_sig, strata$n_total))
strata$any_lo  <- 100 * ci2[,1]; strata$any_hi  <- 100 * ci2[,2]

long <- bind_rows(
  strata %>% transmute(stratum, metric = "Concordant (sig + same direction)",
    pct = pct_concordant, n_pos = n_concordant, n_total, lo = conc_lo, hi = conc_hi),
  strata %>% transmute(stratum, metric = "Any-direction sig overlap",
    pct = pct_our_sig, n_pos = n_our_sig, n_total, lo = any_lo, hi = any_hi)
) %>% mutate(stratum = factor(stratum, levels = c("Unbound","Bound")),
             metric = factor(metric, levels = c("Concordant (sig + same direction)",
                                                 "Any-direction sig overlap")))

# ---- one-sided Fisher's exact test (Bound > Unbound) per metric (identical to source) ----
make_p <- function(metric) {
  m <- strata
  bound   <- m[m$stratum == "Bound", ]
  unbound <- m[m$stratum == "Unbound", ]
  if (metric == "Concordant (sig + same direction)") {
    a <- bound$n_concordant;   b <- bound$n_total - a
    c_ <- unbound$n_concordant; d <- unbound$n_total - c_
  } else {
    a <- bound$n_our_sig; b <- bound$n_total - a
    c_ <- unbound$n_our_sig; d <- unbound$n_total - c_
  }
  ft <- fisher.test(matrix(c(a,b,c_,d), 2, byrow = TRUE), alternative = "greater")
  list(p = ft$p.value, OR = ft$estimate)
}
ann <- data.frame(metric = levels(long$metric), stringsAsFactors = FALSE)
ann$p  <- sapply(ann$metric, function(M) make_p(M)$p)
ann$OR <- sapply(ann$metric, function(M) make_p(M)$OR)
ann$label <- sprintf("OR = %.2f, P = %s", ann$OR,
                     ifelse(ann$p < 0.001, formatC(ann$p, "e", digits=1),
                            sprintf("%.3f", ann$p)))
ann$y_top <- 60

palette <- c("Unbound" = "grey60", "Bound" = "#d73027")

# ---- build the plot (original parameters preserved) ----
pg <- ggplot(long, aes(x = stratum, y = pct, fill = stratum)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct, n_pos, n_total)),
            vjust = -0.3, size = 3, lineheight = 0.85) +
  geom_text(data = ann, aes(x = 1.5, y = y_top, label = label),
            inherit.aes = FALSE, size = 3.2, colour = "grey25") +
  scale_fill_manual(values = palette, guide = "none") +
  facet_grid(metric ~ .) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 60, 20),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "HNRNPM-bound DE genes replicate better between Ho 2021 and our PSD58 limma",
       subtitle = "Stratum = Ho strict-sig restricted to our test universe; error bars = Wilson 95% CI",
       x = NULL, y = "% Ho strict-sig replicating") +
  theme_bw(base_size = 11) +
  theme(plot.title       = element_text(size = 11),                       # original proportion
        plot.subtitle    = element_text(size = 9, colour = "grey40"),     # original proportion
        strip.background = element_rect(fill = "grey90", colour = NA),
        axis.title       = element_text(size = axis_render_pt),   # -> 8 pt after shrink
        axis.text        = element_text(size = axis_render_pt))   # -> 8 pt after shrink

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_Sz.pdf"), pg,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_Sz.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_Sz.pdf)
# ============================================================
# Replication of the Ho et al. 2021 HNRNPM-knockdown DE genes in our additive
# (~ group + PSD) limma analysis, stratified by whether each gene is HNRNPM-bound in
# the Ho CLIP data (red) or unbound (grey). A Ho DE gene is defined by the stringent
# criterion (q < 0.05 in BOTH clones B7 and B9, same direction of change), restricted
# to the genes present in our test universe. The top facet counts a gene as replicated
# only if it is significant in our data AND changes in the same direction as Ho
# (concordant); the bottom facet counts any significant overlap regardless of direction.
# Bars show the percentage of Ho strict-sig genes in each stratum that replicate
# (n replicated / n in stratum printed above each bar), with the percentage denominator
# taken within each stratum; error bars are Wilson 95% confidence intervals. Replication
# of HNRNPM-bound versus unbound genes was compared by a one-sided Fisher's exact test
# (alternative: Bound > Unbound); the odds ratio and P value are annotated on each facet.
# The horizontal axis is the binding stratum and the vertical axis is the percentage of
# Ho strict-sig genes replicating. HNRNPM-bound targets replicate at a modestly but
# significantly higher rate than unbound genes, consistent with direct HNRNPM targets
# being less subject to cell-type-specific rewiring.
