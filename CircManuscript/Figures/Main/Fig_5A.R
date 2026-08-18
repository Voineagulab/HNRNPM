# Fig 5A — SH-SY5Y growth curve, gHNRNPM_1 vs gCtrl.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w2_growth_curve/growth_curve.R
# The source draws two experiments (gHNRNPM_1, gHNRNPM_2) as two PDF pages; this
#   reproduces ONLY the gHNRNPM_1 growth curve (sheet GrowthCurve_gHNRNPM_1_raw),
#   recomputing the mixed two-way RM ANOVA (Greenhouse-Geisser) and per-day
#   Bonferroni t-tests exactly as the source, and drawing the same curve.
#
# Strategy (same as the Fig_2/Fig_3/Fig_4 series): keep ALL original plot parameters
#   (line/point/errorbar, ANOVA annotation text size, significance-star size, title),
#   render at the original size, and let ggsave(scale=) squash the figure to
#   7.75 x 6.83 cm (uniform shrink preserves proportions). Pre-inflate ONLY the axis
#   title/tick fonts + legend so they land at their target pt size in the output
#   (axis titles and legend at 8 pt, x/y tick labels at 4 pt; see TICK_SCALE).
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger
#   than nominal; resize it to 7.75 cm wide on the slide. See Fig_2E_exactsize.R.
# DEVIATION from the source: the y-axis tick labels are drawn in scientific notation
#   with no decimal places, e.g. 2e+05 (sci_labels() below), instead of scales::comma.
#   The x and y tick labels are also set to half the original axis font size, i.e.
#   4 pt final instead of 8 pt. Axis titles and the legend are unchanged at 8 pt.
#   The four-row ANOVA block sits in the bottom-left corner of the panel (it was
#   top-left), and the legend stays to the right of the panel as in the original but
#   is bottom-aligned with the foot of the y axis (it was vertically centred).

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(rstatix); library(scales)
})


# --- Path definition ---
GC_DIR     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w2_growth_curve"
INPUT_XLSX <- file.path(GC_DIR, "growth_curve.xlsx")
SHEET      <- "GrowthCurve_gHNRNPM_1_raw"
LABEL      <- "gHNRNPM_1"
OUT_DIR    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/Figures/Main"

BASE_FONT  <- 11
days       <- 0:6

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 6                     # source page width = FIG_WIDTH_MM/20 = 120/20 = 6 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis title + legend font size
TICK_SCALE <- 0.5                  # x/y tick labels at 50% of AXIS_PT -> 4 pt final
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink
tick_render_pt <- AXIS_PT * TICK_SCALE * SC   # ditto, landing at 4 pt after the shrink

# ggplot geom text sizes are in mm (pt / 2.845); kept at the source's original values
anova.size <- (BASE_FONT - 2) / 2.845
star.size  <- BASE_FONT / 2.845

# ---- y-axis tick labels: scientific notation, no decimal places ----
# 200000 -> "2e+05". sprintf("%.0e") does the mantissa rounding and the exponent
#   bump itself (99999 -> "1e+05"), so no manual normalisation is needed here.
# Exact zero is kept as "0" rather than "0e+00", and NA breaks render blank.
# CAVEAT: with no decimals a tick only reads exactly when the break is a whole
#   multiple of a power of ten. The breaks on this panel are 0, 2e5, 4e5 and 6e5
#   (as in the previous Fig_5A.pdf), so all four labels are exact. If new data
#   ever moves the breaks onto 1.5e5 steps, two ticks would both read "2e+05";
#   in that case add breaks = seq(0, 6e5, by = 2e5) to scale_y_continuous().
sci_labels <- function(x) {
  out <- sprintf("%.0e", x)
  out[is.na(x)] <- ""
  out[!is.na(x) & x == 0] <- "0"
  out
}

# ---- Read the wide sheet, reshape to long ----
wide <- as.data.frame(read_excel(INPUT_XLSX, sheet = SHEET))
long <- data.frame()
for (r in 1:nrow(wide)) for (d in days)
  long <- rbind(long, data.frame(
    subject   = paste0(wide$Cell_line[r], "_", wide$Replicate[r]),
    treatment = wide$Cell_line[r], time = d, count = wide[r, paste0("D", d)]))
long$treatment <- factor(long$treatment, levels = c("gCtrl", "gHNRNPM"))
long$subject   <- factor(long$subject)
long$day       <- factor(long$time)

# ---- Mixed two-way RM ANOVA (Greenhouse-Geisser) ----
fit     <- anova_test(data = long, dv = count, wid = subject, between = treatment, within = day)
tab.gg  <- get_anova_table(fit, correction = "GG")
p.treat <- tab.gg$p[tab.gg$Effect == "treatment"]
p.time  <- tab.gg$p[tab.gg$Effect == "day"]
p.inter <- tab.gg$p[tab.gg$Effect == "treatment:day"]

# ---- Per-day gCtrl vs gHNRNPM (pooled t-test, Bonferroni over Day 1-6) ----
test.days <- 1:6
p.raw <- sapply(test.days, function(d)
  t.test(count ~ treatment, data = long[long$time == d, ], var.equal = TRUE)$p.value)
p.adj <- p.adjust(p.raw, method = "bonferroni")

# ---- mean +/- SD per group/day ----
mean.tab <- aggregate(count ~ treatment + time, long, mean)
sd.tab   <- aggregate(count ~ treatment + time, long, sd)
summ <- mean.tab; colnames(summ)[colnames(summ) == "count"] <- "mean"; summ$sd <- sd.tab$count

# ---- star positions ----
star.df <- data.frame()
for (k in seq_along(test.days)) {
  d <- test.days[k]; p <- p.adj[k]
  star <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
  if (star != "") {
    this.day <- summ[summ$time == d, ]
    star.df <- rbind(star.df, data.frame(time = d, star = star, y = max(this.day$mean + this.day$sd) * 1.05))
  }
}

# ---- build the growth-curve plot (original parameters preserved) ----
y.max <- max(summ$mean + summ$sd)
anova.text <- paste0("Two-way RM ANOVA (Greenhouse-Geisser)\n",
                     "Treatment: p = ", signif(p.treat, 3), "\n",
                     "Time: p = ", signif(p.time, 3), "\n",
                     "Treatment x Time: p = ", signif(p.inter, 3))

p <- ggplot(summ, aes(x = time, y = mean, color = treatment)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
  # ANOVA block in the bottom-left corner of the panel: anchored at Day 0 (hjust = 0)
  # and sitting on the x-axis (vjust = 0 grows the four lines upward from y = 0).
  # y = 0 rather than a small positive offset because the four lines span roughly
  # 110,000 cells on this scale and the D0 markers sit at 120,000; lifting the block
  # even 2% of y.max makes its top line collide with them.
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,
           size = anova.size, label = anova.text) +
  scale_color_manual(values = c("gCtrl" = "grey40", "gHNRNPM" = "orange"),
                     labels = c("gCtrl", LABEL), name = "") +
  scale_x_continuous(breaks = days, labels = paste0("D", days)) +
  scale_y_continuous(labels = sci_labels, limits = c(0, y.max * 1.32)) +
  ggtitle(paste0("Growth curve: ", LABEL, " vs gCtrl")) +
  xlab("Day") + ylab("Cell number") +
  theme_classic(base_size = BASE_FONT) +
  theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
        axis.title   = element_text(size = axis_render_pt),  # -> 8 pt after shrink
        axis.text    = element_text(size = tick_render_pt),  # -> 4 pt after shrink
        legend.text  = element_text(size = axis_render_pt),  # -> 8 pt after shrink
        legend.title = element_text(size = axis_render_pt),
        # Legend outside the panel on its right, as in the original, but bottom-
        # aligned rather than vertically centred. For an outside legend,
        # justification "bottom" anchors the guide box to the bottom of the PANEL,
        # so with a zero bottom margin the legend's lower edge lands level with the
        # foot of the y axis and never dips below it. Give the margin a positive b
        # value to float the legend above the axis instead.
        legend.position      = "right",
        legend.direction     = "vertical",
        legend.justification = "bottom",
        legend.margin        = margin(t = 0, r = 0, b = 0, l = 0))
if (nrow(star.df) > 0)
  p <- p + geom_text(data = star.df, aes(x = time, y = y, label = star),
                     color = "black", size = star.size, inherit.aes = FALSE)

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_5A.pdf"), p,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_5A.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_5A.pdf)
# ============================================================
# Proliferation of SH-SY5Y (dCas9-KRAB) cells over six days after CRISPRi knockdown
# of HNRNPM with guide gHNRNPM_1, compared to a non-targeting control (gCtrl). Points
# are the mean cell number per well and error bars are +/- standard deviation across
# triplicate wells (n = 3); lines connect the daily means for each group (gCtrl in
# grey, gHNRNPM_1 in orange). The horizontal axis is day (D0 seeding to D6) and the
# vertical axis is cell number. The in-panel text reports a mixed two-way
# repeated-measures ANOVA with Greenhouse-Geisser correction (treatment as the
# between-subjects factor, time as the within-subjects factor). Asterisks mark days
# at which gCtrl and gHNRNPM_1 differ by a pooled two-sample t-test after Bonferroni
# correction across days 1 to 6 (* P < 0.05, ** P < 0.01, *** P < 0.001).
