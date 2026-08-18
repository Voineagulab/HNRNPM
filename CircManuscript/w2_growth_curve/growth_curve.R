# growth_curve.R
#
# Author: Juli Wang
#
# Purpose
# -------
# Analyse and plot the SH-SY5Y (dCas9-KRAB) growth-curve cell-count data in
# growth_curve.xlsx. Two knock-down experiments are analysed, one per sheet:
#   gHNRNPM_1 -> sheet "GrowthCurve_gHNRNPM_1_raw"
#   gHNRNPM_2 -> sheet "GrowthCurve_gHNRNPM_2_raw"
# Each sheet holds triplicate counts (rep1-3) for a non-targeting control
# (gCtrl) and an HNRNPM-targeting guide (gHNRNPM) over Day 0 to Day 6.
#
# For each experiment the script:
#   1. Runs a mixed two-way repeated-measures ANOVA: treatment (gCtrl vs
#      gHNRNPM) is the between-subjects factor and time (Day 0-6) the
#      within-subjects factor, with the Greenhouse-Geisser correction.
#   2. Compares gCtrl vs gHNRNPM at each day (Day 1-6) with a pooled two-sample
#      t-test, Bonferroni-corrected across the 6 days.
#   3. Draws a growth curve (mean +/- SD, n = 3) with the ANOVA result written
#      in the panel and significance stars above the significant days.
#
# Outputs
# -------
#   OUTPUT_PDF : the two growth-curve plots, one per page.
#   OUTPUT_CSV : every statistic computed, in one tidy table, so that reading
#                this script together with the CSV tells you exactly what was
#                done and what to report - the descriptive mean/SD per group and
#                day, the per-day t-tests (raw and Bonferroni p), and the overall
#                ANOVA effects (F, df, Greenhouse-Geisser epsilon, effect size,
#                p). Each row also has a ready-to-quote "report" string.
#
# Run this in RStudio (Source), or with: Rscript growth_curve.R

## ----------------------------- settings ---------------------------------- ##
INPUT_XLSX  <- "growth_curve.xlsx"
OUTPUT_PDF  <- "growth_curve.pdf"
OUTPUT_CSV  <- "stats.csv"
RESULTS_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w2_growth_curve/"
setwd(RESULTS_DIR)

# One experiment per plot: sheet name -> plot label
EXPERIMENTS <- list(
  "GrowthCurve_gHNRNPM_1_raw" = "gHNRNPM_1",
  "GrowthCurve_gHNRNPM_2_raw" = "gHNRNPM_2"
)

# Manuscript figure size (mm) per page; converted to inches for pdf().
FIG_WIDTH_MM  <- 120
FIG_HEIGHT_MM <- 95
BASE_FONT     <- 11     # base font size (pt)

days <- 0:6            # the time points, Day 0 to Day 6

## ----------------------------- packages ---------------------------------- ##
library(readxl)     # read_excel(): read the .xlsx sheets
library(ggplot2)    # plotting
library(rstatix)    # anova_test(): mixed RM ANOVA with Greenhouse-Geisser

# ggplot2 measures text "size" in mm, so divide a point size by 2.845 to convert.
anova.size <- (BASE_FONT - 2) / 2.845   # slightly smaller text for the ANOVA box
star.size  <- BASE_FONT / 2.845         # significance stars

## ---------------- loop over the two experiments --------------------------- ##
plot.list <- list()      # holds the two growth-curve plots
stats.all <- data.frame()  # collects every statistic; written to OUTPUT_CSV at the end

for (i in 1:length(EXPERIMENTS)) {
  sheet <- names(EXPERIMENTS)[i]
  label <- EXPERIMENTS[[i]]

  ## --- read the wide sheet and reshape to long (one row per replicate/day) --
  wide <- as.data.frame(read_excel(INPUT_XLSX, sheet = sheet))

  long <- data.frame()
  for (r in 1:nrow(wide)) {
    for (d in days) {
      one.row <- data.frame(
        subject   = paste0(wide$Cell_line[r], "_", wide$Replicate[r]),
        treatment = wide$Cell_line[r],
        time      = d,
        count     = wide[r, paste0("D", d)])
      long <- rbind(long, one.row)
    }
  }

  long$treatment <- factor(long$treatment, levels = c("gCtrl", "gHNRNPM"))
  long$subject   <- factor(long$subject)
  long$day       <- factor(long$time)   # factor version of time, needed by the ANOVA

  ## --- mixed two-way RM ANOVA (treatment between, time within) --------------
  # Get the table twice: uncorrected and Greenhouse-Geisser corrected. Dividing
  # the GG numerator df by the uncorrected numerator df gives the GG epsilon.
  fit      <- anova_test(data = long, dv = count, wid = subject,
                         between = treatment, within = day)
  tab.none <- get_anova_table(fit, correction = "none")
  tab.gg   <- get_anova_table(fit, correction = "GG")

  p.treat <- tab.gg$p[tab.gg$Effect == "treatment"]
  p.time  <- tab.gg$p[tab.gg$Effect == "day"]
  p.inter <- tab.gg$p[tab.gg$Effect == "treatment:day"]

  ## --- post-hoc: gCtrl vs gHNRNPM at each day (Day 1-6) ---------------------
  # Pooled two-sample t-test on each day, then Bonferroni across the 6 days.
  # Day 0 is the seeding day and is not compared.
  test.days <- 1:6
  p.raw <- c()
  for (d in test.days) {
    sub <- long[long$time == d, ]
    p.raw <- c(p.raw, t.test(count ~ treatment, data = sub, var.equal = TRUE)$p.value)
  }
  p.adj <- p.adjust(p.raw, method = "bonferroni")

  ## --- print the numbers to the console ------------------------------------
  print(paste0("=============== ", label, " vs gCtrl ==============="))
  print(paste0("Treatment p = ", signif(p.treat, 3)))
  print(paste0("Time p = ", signif(p.time, 3)))
  print(paste0("Treatment x Time p = ", signif(p.inter, 3)))
  print(data.frame(day = test.days, p_bonferroni = signif(p.adj, 3)))

  ## --- collect the ANOVA effects into the stats table ----------------------
  # Nice names for the three effects (keyed by the ANOVA table's Effect names).
  term.name <- c("treatment"     = "Treatment (between-subjects)",
                 "day"           = "Time (within-subjects)",
                 "treatment:day" = "Treatment x Time (within-subjects)")

  for (j in 1:nrow(tab.gg)) {
    eff  <- tab.gg$Effect[j]
    Fval <- tab.gg$F[j]
    df1  <- tab.gg$DFn[j]                      # GG-corrected numerator df
    df2  <- tab.gg$DFd[j]                      # GG-corrected denominator df
    pval <- tab.gg$p[j]
    ges  <- tab.gg$ges[j]                      # generalized eta-squared (effect size)
    eps  <- tab.gg$DFn[j] / tab.none$DFn[j]    # Greenhouse-Geisser epsilon (1 = between effect)
    sig  <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns")))

    report <- paste0("F(", round(df1, 2), ",", round(df2, 2), ") = ", signif(Fval, 3),
                     ", p = ", signif(pval, 3), " (GG epsilon = ", round(eps, 3), ")")

    stats.all <- rbind(stats.all, data.frame(
      experiment = label,
      analysis   = "Mixed two-way RM ANOVA (Greenhouse-Geisser)",
      term       = term.name[eff],
      mean_gCtrl = NA, sd_gCtrl = NA, n_gCtrl = NA,
      mean_gHNRNPM = NA, sd_gHNRNPM = NA, n_gHNRNPM = NA,
      statistic  = signif(Fval, 4),
      df1 = round(df1, 3), df2 = round(df2, 3),
      GG_epsilon = round(eps, 3),
      effect_size_ges = signif(ges, 3),
      p_value = signif(pval, 3),
      p_adjust_method = NA, p_adjusted = NA,
      significance = sig,
      report = report))
  }

  ## --- collect the per-day descriptives + t-tests into the stats table -----
  for (d in days) {
    ct <- long$count[long$time == d & long$treatment == "gCtrl"]
    hn <- long$count[long$time == d & long$treatment == "gHNRNPM"]

    if (d %in% test.days) {
      tt   <- t.test(count ~ treatment, data = long[long$time == d, ], var.equal = TRUE)
      praw <- p.raw[d]           # test.days is 1:6, so day d sits at index d
      padj <- p.adj[d]
      sig  <- ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "ns")))
      report <- paste0("t(", round(unname(tt$parameter), 1), ") = ", signif(unname(tt$statistic), 3),
                       ", p_raw = ", signif(praw, 3), ", p_bonferroni = ", signif(padj, 3))

      stats.all <- rbind(stats.all, data.frame(
        experiment = label,
        analysis   = "Per-day gCtrl vs gHNRNPM (pooled two-sample t-test)",
        term       = paste0("Day ", d),
        mean_gCtrl = round(mean(ct)), sd_gCtrl = round(sd(ct)), n_gCtrl = length(ct),
        mean_gHNRNPM = round(mean(hn)), sd_gHNRNPM = round(sd(hn)), n_gHNRNPM = length(hn),
        statistic  = signif(unname(tt$statistic), 4),
        df1 = unname(tt$parameter), df2 = NA,
        GG_epsilon = NA, effect_size_ges = NA,
        p_value = signif(praw, 3),
        p_adjust_method = "bonferroni", p_adjusted = signif(padj, 3),
        significance = sig,
        report = report))
    } else {
      # Day 0: seeding day, described but not tested.
      stats.all <- rbind(stats.all, data.frame(
        experiment = label,
        analysis   = "Descriptive only (seeding day, not compared)",
        term       = paste0("Day ", d),
        mean_gCtrl = round(mean(ct)), sd_gCtrl = round(sd(ct)), n_gCtrl = length(ct),
        mean_gHNRNPM = round(mean(hn)), sd_gHNRNPM = round(sd(hn)), n_gHNRNPM = length(hn),
        statistic  = NA, df1 = NA, df2 = NA,
        GG_epsilon = NA, effect_size_ges = NA,
        p_value = NA, p_adjust_method = NA, p_adjusted = NA,
        significance = NA,
        report = "Day 0 seeding; identical plating, not tested"))
    }
  }

  ## --- mean and SD per group/day for the plot ------------------------------
  mean.tab <- aggregate(count ~ treatment + time, long, mean)
  sd.tab   <- aggregate(count ~ treatment + time, long, sd)
  summ <- mean.tab
  colnames(summ)[colnames(summ) == "count"] <- "mean"
  summ$sd <- sd.tab$count

  ## --- work out where to put the significance stars ------------------------
  # A star sits just above the taller (mean + SD) of the two groups that day.
  star.df <- data.frame()
  for (k in 1:length(test.days)) {
    d <- test.days[k]
    p <- p.adj[k]

    star <- ""
    if (p < 0.001) { star <- "***"
    } else if (p < 0.01) { star <- "**"
    } else if (p < 0.05) { star <- "*" }

    if (star != "") {
      this.day <- summ[summ$time == d, ]
      y.top <- max(this.day$mean + this.day$sd)
      star.df <- rbind(star.df, data.frame(time = d, star = star, y = y.top * 1.05))
    }
  }

  ## --- build the growth-curve plot -----------------------------------------
  y.max <- max(summ$mean + summ$sd)

  anova.text <- paste0("Two-way RM ANOVA (Greenhouse-Geisser)\n",
                       "Treatment: p = ", signif(p.treat, 3), "\n",
                       "Time: p = ", signif(p.time, 3), "\n",
                       "Treatment x Time: p = ", signif(p.inter, 3))

  p <- ggplot(summ, aes(x = time, y = mean, color = treatment)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
    annotate("text", x = 0, y = y.max * 1.28, hjust = 0, vjust = 1,
             size = anova.size, label = anova.text) +
    scale_color_manual(values = c("gCtrl" = "grey40", "gHNRNPM" = "orange"),
                       labels = c("gCtrl", label), name = "") +
    scale_x_continuous(breaks = days, labels = paste0("D", days)) +
    scale_y_continuous(labels = scales::comma, limits = c(0, y.max * 1.32)) +
    ggtitle(paste0("Growth curve: ", label, " vs gCtrl")) +
    xlab("Day") + ylab("Cell number") +
    theme_classic(base_size = BASE_FONT) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # add the stars (inherit.aes = FALSE so they are not coloured by treatment)
  if (nrow(star.df) > 0) {
    p <- p + geom_text(data = star.df, aes(x = time, y = y, label = star),
                       color = "black", size = star.size, inherit.aes = FALSE)
  }

  plot.list[[i]] <- p
}

## ------------------------- save the outputs ------------------------------- ##
# Statistics table: one row per descriptive/per-day test/ANOVA effect.
write.csv(stats.all, OUTPUT_CSV, row.names = FALSE)
print(paste0("Saved ", OUTPUT_CSV))

# Plots: one experiment per page, both in the single output PDF.
pdf(OUTPUT_PDF, width = FIG_WIDTH_MM / 20, height = FIG_HEIGHT_MM / 25.4)
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()
print(paste0("Saved ", OUTPUT_PDF))
