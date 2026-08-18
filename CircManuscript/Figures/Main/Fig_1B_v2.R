# Fig 1B v2 — HnrnpM eCLIP enrichment at circRNA flanking introns (box plot),
#             with four flank-vs-background comparisons drawn as long dashes.
#
# Exports ONE file: Figures/Main/Fig_1B_v2.pdf
#
# Base plot adapted from: Figures/Main/Fig_1B.R (itself panel A of
#   ForPublication/FinalAnalysis/SCRIPTS/09_HnrnpM_CLIP_enrichment/plot_enrichment_PSD58.R).
# Boxes 1-4 and 6, the factor order, palette, theme, fonts and physical size are
# UNCHANGED from Fig_1B.R. Three things are added or changed:
#   1. box 5 is REDRAWN (see below),
#   2. a fixed -10 to 10 y axis (which clips nothing — see below),
#   3. four solid horizontal comparison bars, each with a one-sided P value above it.
#
# ---------------------------------------------------------------------------
# TERMINOLOGY — what a "region" is
# One row of clip_signal_per_region.tsv = one genomic interval (an intron) with its own
# IP and Input coverage, giving one log2(IP/Input) value = one point in a box.
#   flank5          exactly 1 interval per circRNA (the intron 5' of the BSJ)
#   flank3          exactly 1 interval per circRNA (the intron 3' of the BSJ)
#   within_gene_bg  5 intervals per circRNA (other, length-matched introns of the SAME
#                   host gene). 55 DE circRNAs have these: 54 x 5 + 1 x 3 = 273 intervals.
#   non_host_bg     4,999 intervals sampled from expressed genes that host no circRNA;
#                   not tied to any circRNA, so it can never be paired.
#
# ---------------------------------------------------------------------------
# BOX 5 IS NOT THE BOX Fig_1B.R DRAWS
# Fig_1B.R labels box 5 "Within-gene bg" and fills it with the within-gene background of
# the WHOLE circRNA universe — 27,614 intervals from 5,570 circRNAs, of which only 273
# (1.0%, 55 circRNAs) come from the DE set. Its de_lab prefix is applied to flank5/flank3
# only, never to within_gene_bg, so that box is ~99% non-DE. Comparing DE flanks against
# it measures partly a gene-set difference: DE host genes carry more HNRNPM signal than
# average genes (their own background sits at -0.24 vs -0.59 for the pooled box).
#
# Here box 5 is rebuilt as the within-gene background OF THE DE circRNAs ONLY, collapsed
# to one value per circRNA by averaging that circRNA's 5 background introns — 55 points,
# relabelled "DE within-gene bg". Collapsing is unavoidable for a paired test: pairing
# needs exactly one background number per circRNA. Each DE circRNA is then its own
# control, so host-gene binding level cancels out and the flank-specific effect is
# isolated. The resulting box median (-0.2394) is exactly the median2 that
# enrichment_test_summary.tsv already publishes for BSJ_additive test 2.
#
# CAVEAT to state in the legend: averaging 5 intervals shrinks this box's spread
# (SD 1.12) relative to the single-interval flank boxes (SD 1.52 and 1.32). The tighter
# box is an artefact of averaging, not biology.
#
# ---------------------------------------------------------------------------
# COMPARISONS DRAWN (bottom bar to top bar)
#   DE 3' flank vs DE within-gene bg     (box 2 vs box 5)   paired
#   DE 5' flank vs DE within-gene bg     (box 1 vs box 5)   paired
#   DE 3' flank vs Non-host-gene bg      (box 2 vs box 6)   unpaired
#   DE 5' flank vs Non-host-gene bg      (box 1 vs box 6)   unpaired
# Ordered shortest span at the bottom so the bars nest without colliding.
#
# ---------------------------------------------------------------------------
# STATISTICAL TEST — matched to enrichment_test_summary.tsv
# That table was produced by test_CLIP_enrichment_PSD58.py using scipy:
#   * scipy.stats.wilcoxon(alternative="greater")     paired, its test 2 (within-gene bg)
#   * scipy.stats.mannwhitneyu(alternative="greater") unpaired, its test 3 (non-host bg)
# At these n scipy resolves method="auto" to the ASYMPTOTIC (normal-approximation) form
# in every case (verified against method="asymptotic"). The R equivalents are pinned:
#   paired   : wilcox.test(paired = TRUE, exact = FALSE, correct = FALSE)
#                                                           -> scipy correction=False
#   unpaired : wilcox.test(exact = FALSE, correct = TRUE)   -> scipy use_continuity=True
# Both drop zero differences (R and scipy zero_method="wilcox" agree; there is exactly
# one zero difference in each paired set). P values are ONE-SIDED (greater), matching the
# p_one_sided column and the Results text, and RAW — the summary table is uncorrected too.
# Verified: all four P values reproduce the original scipy run to 4 significant figures.
#
# The box shows all 55 DE circRNAs that have a within-gene background; the paired tests
# use the 47 / 48 of those that ALSO have the relevant flanking intron.
# ---------------------------------------------------------------------------

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES09   <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/09_HnrnpM_CLIP_enrichment")
SIG     <- file.path(RES09, "clip_signal_per_region.tsv")
BSJ_ADD <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")
OUT_PDF <- file.path(OUT_DIR, "Fig_1B_v2.pdf")      # the ONLY file this script writes

# ---- Target physical size on the slide (unchanged from Fig_1B.R) ----
W_CM   <- 5.77; H_CM <- 6.83      # width x height as it should appear in PowerPoint
MIN_PT <- 8                        # font size (points)

# ---- Axis limits and bracket geometry (data units on the y scale) ----
# Fixed -10 to 10 axis. With box 5 restricted to the DE circRNAs the data span only
# -6.31 to +6.09 (the +9.02 outliers lived in the pooled universe-wide background that
# this figure no longer draws), so these limits CLIP NOTHING — every point, whisker and
# outlier is inside the frame.
Y_LO     <- -10                    # axis floor  (data min -6.31)
Y_HI     <-  10                    # axis ceiling (data max +6.09)
# Comparison bars sit above the tallest whisker (+3.70). Step 1.5 leaves room for an
# 8 pt P value drawn ABOVE each bar without it touching the bar above; the topmost label
# lands near +9.8, inside the frame.
BR_Y0    <-  3.9                   # first (lowest) bar
BR_STEP  <-  1.5                   # vertical spacing between bars
LAB_MM   <- MIN_PT / .pt           # geom_text size in mm that renders at MIN_PT points

BG_LAB   <- "DE within-gene bg"    # box 5 label (renamed: no longer the universe-wide bg)

# ---- Load data ----
sig <- read.table(SIG,     sep = "\t", header = TRUE, check.names = FALSE)
bsj <- read.table(BSJ_ADD, sep = "\t", header = TRUE, check.names = FALSE)

# Additive BSJ-DE circRNAs (publication primary)
de_circs <- bsj$circRNA[bsj$adj.P.Val < 0.05]

# ---- Panel labels ----
sig$de_lab <- ifelse(sig$circRNA %in% de_circs, "DE", "non-DE")
sig$panel_label <- with(sig, ifelse(
  region_type == "flank5", paste0(de_lab, " 5' flank"),
  ifelse(region_type == "flank3", paste0(de_lab, " 3' flank"),
  ifelse(region_type == "non_host_bg", "Non-host-gene bg", NA_character_))))

LEVELS  <- c("DE 5' flank","DE 3' flank","non-DE 5' flank","non-DE 3' flank",
             BG_LAB,"Non-host-gene bg")
palette <- c("DE 5' flank" = "#d73027", "DE 3' flank" = "#d73027",
             "non-DE 5' flank" = "#4C78A8", "non-DE 3' flank" = "#4C78A8",
             "Non-host-gene bg" = "grey75")
palette[BG_LAB] <- "grey50"        # same grey as the box it replaces

# ---- Box 5: DE-only within-gene background, one averaged value per circRNA ----
wg_de <- subset(sig, region_type == "within_gene_bg" &
                     circRNA %in% de_circs & !is.na(log2_IP_over_Input))
wg_de_percirc <- aggregate(log2_IP_over_Input ~ circRNA, data = wg_de, FUN = mean)
names(wg_de_percirc)[2] <- "wg_score"

# ---- Assemble the plotting frame (boxes 1-4 and 6 verbatim, box 5 rebuilt) ----
oth <- subset(sig, !is.na(panel_label) & !is.na(log2_IP_over_Input))
dat <- rbind(
  data.frame(panel_label = oth$panel_label,     value = oth$log2_IP_over_Input),
  data.frame(panel_label = BG_LAB,              value = wg_de_percirc$wg_score)
)
dat$panel_label <- factor(dat$panel_label, levels = LEVELS)

# ===========================================================================
# STATISTICS
# ===========================================================================

de5    <- sig$log2_IP_over_Input[sig$region_type == "flank5" & sig$circRNA %in% de_circs]
de3    <- sig$log2_IP_over_Input[sig$region_type == "flank3" & sig$circRNA %in% de_circs]
nh_all <- sig$log2_IP_over_Input[sig$region_type == "non_host_bg"]
de5 <- de5[!is.na(de5)]; de3 <- de3[!is.na(de3)]; nh_all <- nh_all[!is.na(nh_all)]

# Pair each DE circRNA's flanking intron with its own averaged within-gene background
paired_vals <- function(rt) {
  f <- subset(sig, region_type == rt & circRNA %in% de_circs & !is.na(log2_IP_over_Input),
              select = c("circRNA", "log2_IP_over_Input"))
  m <- merge(f, wg_de_percirc, by = "circRNA")
  m[!is.na(m$wg_score), ]
}
pair5 <- paired_vals("flank5")
pair3 <- paired_vals("flank3")

# scipy-matched wrappers (see header note)
p_unpaired <- function(x, y)
  wilcox.test(x, y, alternative = "greater", exact = FALSE, correct = TRUE)$p.value
p_paired <- function(x, y)
  wilcox.test(x, y, paired = TRUE, alternative = "greater", exact = FALSE, correct = FALSE)$p.value

p_5_wg <- p_paired(pair5$log2_IP_over_Input, pair5$wg_score)
p_3_wg <- p_paired(pair3$log2_IP_over_Input, pair3$wg_score)
p_5_nh <- p_unpaired(de5, nh_all)
p_3_nh <- p_unpaired(de3, nh_all)

# Pretty p-value formatter (same convention as Fig_2Z.R)
fp <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "g", digits = 2)
}

# ===========================================================================
# PLOT
# ===========================================================================

br <- data.frame(
  x    = c(2, 1, 2, 1),                        # DE 3'|bg, DE 5'|bg, DE 3'|nh, DE 5'|nh
                                               # (shortest span first = lowest bar)
  xend = c(5, 5, 6, 6),
  y    = BR_Y0 + BR_STEP * (0:3),
  lab  = paste0("P = ", vapply(c(p_3_wg, p_5_wg, p_3_nh, p_5_nh), fp, character(1))),
  stringsAsFactors = FALSE
)

p <- ggplot(dat, aes(x = panel_label, y = value, fill = panel_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, colour = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = palette, guide = "none") +
  geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linetype = "solid",
               linewidth = 0.3, colour = "grey20") +
  geom_text(data = br, aes(x = (x + xend) / 2, y = y, label = lab),
            inherit.aes = FALSE, size = LAB_MM, colour = "grey20",
            vjust = -0.4) +                    # P value sits ABOVE its bar
  coord_cartesian(ylim = c(Y_LO, Y_HI)) +
  labs(title = "log2(IP/Input) per region",
       x = NULL, y = "log2(IP / Input)") +
  theme_bw(base_size = MIN_PT) +
  theme(
    text        = element_text(size = MIN_PT),
    axis.text   = element_text(size = MIN_PT),
    axis.title  = element_text(size = MIN_PT),
    plot.title  = element_text(size = MIN_PT, hjust = 0.5, margin = margin(b = 1)),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(OUT_PDF, p, width = W_CM, height = H_CM, units = "cm", device = cairo_pdf)
cat("Wrote:", OUT_PDF, "\n")

# ===========================================================================
# CONSOLE REPORT
# ===========================================================================
cat("\n=== Box descriptives ===\n")
out <- do.call(rbind, lapply(levels(dat$panel_label), function(l) {
  v <- dat$value[dat$panel_label == l]; v <- v[!is.na(v)]
  q <- quantile(v, c(.25, .5, .75))
  data.frame(box = which(levels(dat$panel_label) == l), label = l, n = length(v),
             Q1 = round(q[[1]], 4), median = round(q[[2]], 4), Q3 = round(q[[3]], 4),
             SD = round(sd(v), 3))
}))
print(out, row.names = FALSE)

cat("\n=== Box 5 construction ===\n")
cat(sprintf("  %d intervals from %d DE circRNAs, averaged to %d values (one per circRNA)\n",
            nrow(wg_de), length(unique(wg_de$circRNA)), nrow(wg_de_percirc)))
cat(sprintf("  intervals per DE circRNA: %s\n",
            paste(names(table(table(wg_de$circRNA))), "x",
                  table(table(wg_de$circRNA)), collapse = ", ")))

cat("\n=== Tests (one-sided 'greater', raw P) ===\n")
tt <- data.frame(
  comparison = c("DE 5' flank vs own within-gene bg",
                 "DE 3' flank vs own within-gene bg",
                 "DE 5' flank vs Non-host-gene bg",
                 "DE 3' flank vs Non-host-gene bg"),
  test = c(rep("signed-rank (paired)", 2), rep("rank-sum (unpaired)", 2)),
  n1   = c(nrow(pair5), nrow(pair3), length(de5), length(de3)),
  n2   = c(nrow(pair5), nrow(pair3), length(nh_all), length(nh_all)),
  median1 = round(c(median(pair5$log2_IP_over_Input), median(pair3$log2_IP_over_Input),
                    median(de5), median(de3)), 4),
  median2 = round(c(median(pair5$wg_score), median(pair3$wg_score),
                    median(nh_all), median(nh_all)), 4),
  p_one_sided = signif(c(p_5_wg, p_3_wg, p_5_nh, p_3_nh), 4),
  stringsAsFactors = FALSE
)
print(tt, row.names = FALSE)

cat("\n=== Reference values from the original scipy run (must match above) ===\n")
cat("  DE 5' vs own wg (paired) : 9.459e-06\n")
cat("  DE 3' vs own wg (paired) : 2.302e-06\n")
cat("  DE 5' vs non-host bg     : 7.947e-06\n")
cat("  DE 3' vs non-host bg     : 7.615e-10\n")
cat(sprintf("\n  enrichment_test_summary.tsv BSJ_additive test 2 median2 = -0.2394; box 5 median = %.4f\n",
            median(wg_de_percirc$wg_score)))

cat("\n"); print(sessionInfo()$R.version$version.string)


# ---- Figure legend ----
# (B) Box plot of HNRNPM eCLIP enrichment across six genomic region categories. eCLIP signal is expressed as log2 of immunoprecipitation coverage over input coverage (IP/Input), computed per region from HNRNPM eCLIP data published by Ho et al (2021) lifted to the hg38 assembly. Each region is one intron, contributing one log2(IP/Input) value. The categories are the 5 prime and 3 prime introns flanking the 56 DE circRNAs highlighted in (A), the 5 prime and 3 prime introns flanking non-DE circRNAs, an intronic background drawn from within the host genes of the DE circRNAs, and a background drawn from introns of expressed genes that do not host circRNAs. For the within-host-gene background, five length-matched introns were sampled from each DE circRNA host gene and averaged to give one value per circRNA (55 circRNAs); note that this averaging reduces the spread of that box relative to the single-intron flank boxes. Boxes show the median and interquartile range with whiskers, and the horizontal dashed line marks equal immunoprecipitation and input signal at a log2 value of zero. Higher values indicate stronger HNRNPM binding. Horizontal bars mark four comparisons of the DE flanking introns against the two backgrounds, each annotated above the bar with a one-sided (greater) Wilcoxon P value, uncorrected. Comparisons against the within-host-gene background are paired Wilcoxon signed-rank tests of each DE circRNA's flanking intron against the averaged background of that same host gene (n = 47 for the 5 prime flank and 48 for the 3 prime flank, being the DE circRNAs for which both the relevant flanking intron and a within-gene background were available), so that each circRNA serves as its own control and host-gene binding level cancels out. Comparisons against the non-host-gene background are unpaired Wilcoxon rank-sum tests, since those regions are not associated with any circRNA.
