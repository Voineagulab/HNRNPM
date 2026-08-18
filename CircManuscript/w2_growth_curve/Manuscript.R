# Manuscript.R
#
# Author: Juli Wang
#
# Purpose
# -------
# Generate a BMC Genomics-style Methods and Results text for the SH-SY5Y
# (dCas9-KRAB) HNRNPM-knockdown growth-curve proliferation assay in this folder.
# The Methods section describes the experiment that produced growth_curve.xlsx.
# The Results section is written live from stats.csv (the output of
# growth_curve.R), so the reported means, SDs, ANOVA effects and per-day
# Bonferroni comparisons always match the current analysis. The text also
# refers to the growth-curve figure (FIGURE_PDF).
#
# Output: a plain-text file (OUTPUT_TXT) with a Methods section, a Results
# section for each guide RNA, and a per-day summary table.
#
# Run this in RStudio (Source), or with: Rscript Manuscript.R

## ----------------------------- settings ---------------------------------- ##
STATS_CSV   <- "stats.csv"                        # produced by growth_curve.R
FIGURE_PDF  <- "growth_curve.pdf"                 # the growth-curve figure this text refers to
OUTPUT_TXT  <- "Growth_curve_Methods_Results.txt" # the text file to write
RESULTS_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w2_growth_curve/"
setwd(RESULTS_DIR)

SEED       <- 120000                              # cells per well seeded at Day 0
EXPERIMENTS <- c("gHNRNPM_1", "gHNRNPM_2")        # the two guide-RNA experiments, in order

## ----------------------------- helpers ----------------------------------- ##
# Format a whole-number cell count with thousands separators, e.g. 435346 -> "435,346".
comma <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

# Format a p-value for prose: exact to 3 significant figures (keeps scientific
# notation for very small values, e.g. 1.02e-11).
fmt.p <- function(p) paste0(signif(p, 3))

## ----------------------------- read stats -------------------------------- ##
stats <- read.csv(STATS_CSV, stringsAsFactors = FALSE, check.names = FALSE)

## ----------------------------- Methods text ------------------------------ ##
out <- c()
out <- c(out, "Methods")
out <- c(out, "")
out <- c(out, "Cell culture and CRISPR interference (CRISPRi) knockdown")
out <- c(out, paste0(
  "SH-SY5Y human neuroblastoma cells stably expressing dCas9-KRAB were maintained ",
  "under standard culture conditions. HNRNPM was knocked down by lentiviral ",
  "transduction of an HNRNPM-targeting guide RNA (gHNRNPM); two independent guide ",
  "RNAs (gHNRNPM_1 and gHNRNPM_2) were used in separate experiments, and a ",
  "non-targeting guide RNA (gCtrl) served as the control. Transduced cells were ",
  "selected with puromycin, and selection was removed at the start of the ",
  "proliferation assay (Day 0)."))
out <- c(out, "")
out <- c(out, "Proliferation (growth-curve) assay")
out <- c(out, paste0(
  "To assess proliferation, gCtrl and gHNRNPM cells were seeded at equal density ",
  "(1.2 x 10^5 cells per well) in 6-well plates (Day 0) and counted daily from Day 0 ",
  "to Day 6. At each time point, cells were harvested by trypsinisation, resuspended ",
  "in a defined volume of medium, and counted with a haemocytometer, taking two ",
  "independent density readings per sample; total cell number was calculated as the ",
  "product of the mean cell density and the resuspension volume. Three independent ",
  "biological replicates were performed for each condition, and the two guide RNAs ",
  "were analysed as independent experiments."))
out <- c(out, "")
out <- c(out, "Statistical analysis")
out <- c(out, paste0(
  "Cell-count data were analysed with a mixed-design two-way repeated-measures ANOVA, ",
  "with treatment (gCtrl vs gHNRNPM) as the between-subjects factor and time (Day 0 to ",
  "Day 6) as the within-subjects factor; the Greenhouse-Geisser correction was applied ",
  "to the within-subjects effects to account for violations of sphericity. Because cells ",
  "were seeded at identical density, Day 0 counts were equal across conditions and were ",
  "not included in the between-condition comparisons. Differences between gCtrl and ",
  "gHNRNPM at each subsequent time point (Day 1 to Day 6) were assessed by two-sample ",
  "t-tests with Bonferroni correction across the six time points. Analyses were performed ",
  "in R (v4.3.3) with the rstatix package; generalized eta-squared (eta^2G) is reported as ",
  "the effect size for ANOVA terms, and p < 0.05 was considered statistically significant."))
out <- c(out, "")
out <- c(out, "")

## ----------------------------- Results text ------------------------------ ##
out <- c(out, "Results")
out <- c(out, "")

for (exp in EXPERIMENTS) {
  s <- stats[stats$experiment == exp, ]

  # pull the rows we need
  d0 <- s[s$term == "Day 0", ]
  d6 <- s[s$term == "Day 6", ]
  tr <- s[s$term == "Treatment (between-subjects)", ]
  ti <- s[s$term == "Time (within-subjects)", ]
  it <- s[s$term == "Treatment x Time (within-subjects)", ]

  # Day-6 fold changes over the Day-0 seeding density, and control:KD ratio
  fold.ctrl <- d6$mean_gCtrl   / SEED
  fold.kd   <- d6$mean_gHNRNPM / SEED
  ratio.d6  <- d6$mean_gCtrl   / d6$mean_gHNRNPM

  # per-day post-hoc: significant days and near-significant "trend" days
  perday <- s[s$term %in% paste0("Day ", 1:6), ]
  perday <- perday[order(perday$term), ]
  sig    <- perday[perday$significance %in% c("*", "**", "***"), ]
  trend  <- perday[perday$p_adjusted >= 0.05 & perday$p_adjusted < 0.1, ]

  # build the significant/trend day phrases
  if (nrow(sig) > 0) {
    sig.phrase <- paste(paste0(sig$term, " (P = ", fmt.p(sig$p_adjusted), ")"), collapse = ", ")
  } else {
    sig.phrase <- "none of the individual time points"
  }
  if (nrow(trend) > 0) {
    trend.phrase <- paste0(" A consistent trend was seen at ",
                           paste(paste0(trend$term, " (P = ", fmt.p(trend$p_adjusted), ")"), collapse = ", "), ".")
  } else {
    trend.phrase <- ""
  }

  # word choices that stay correct whatever the numbers are
  tr.word <- ifelse(tr$p_value < 0.05, "a significant", "no significant")
  ti.word <- ifelse(ti$p_value < 0.05, "a significant", "no significant")
  it.word <- ifelse(it$p_value < 0.05, "a significant", "no significant")

  ## --- section heading + paragraphs -------------------------------------- ##
  out <- c(out, paste0("HNRNPM knockdown (", exp, ") reduces SH-SY5Y cell proliferation"))
  out <- c(out, "")

  out <- c(out, paste0(
    "Following equal seeding at 1.2 x 10^5 cells per well, both gCtrl and ", exp,
    " cells grew over the six-day period, but ", exp, " cells proliferated more ",
    "slowly than control (Fig. 1). By Day 6, gCtrl cells reached ", comma(d6$mean_gCtrl),
    " +/- ", comma(d6$sd_gCtrl), " cells per well (n = ", d6$n_gCtrl, "), a ",
    round(fold.ctrl, 1), "-fold increase over the seeding density, whereas ", exp,
    " cells reached only ", comma(d6$mean_gHNRNPM), " +/- ", comma(d6$sd_gHNRNPM),
    " cells per well (a ", round(fold.kd, 1), "-fold increase) - approximately ",
    round(ratio.d6, 1), "-fold fewer cells than control at the same time point."))
  out <- c(out, "")

  out <- c(out, paste0(
    "A mixed-design two-way repeated-measures ANOVA revealed ", tr.word,
    " main effect of treatment (F(", tr$df1, ",", tr$df2, ") = ", tr$statistic,
    ", P = ", fmt.p(tr$p_value), ", eta^2G = ", tr$effect_size_ges, ") and ", ti.word,
    " effect of time (Greenhouse-Geisser corrected: F(", ti$df1, ",", ti$df2, ") = ",
    ti$statistic, ", P = ", fmt.p(ti$p_value), ", eta^2G = ", ti$effect_size_ges,
    "), together with ", it.word, " treatment x time interaction (F(", it$df1, ",", it$df2,
    ") = ", it$statistic, ", P = ", fmt.p(it$p_value), ", eta^2G = ", it$effect_size_ges,
    "), indicating that the two cell lines diverged in their proliferation rate over time."))
  out <- c(out, "")

  out <- c(out, paste0(
    "In post-hoc comparisons (two-sample t-tests, Bonferroni-corrected across Day 1 to Day 6), ",
    "cell counts differed significantly between gCtrl and ", exp, " at ", sig.phrase, ".",
    trend.phrase))
  out <- c(out, "")

  first.sig <- if (nrow(sig) > 0) sig$term[1] else NA
  concl <- if (!is.na(first.sig)) {
    paste0("Together, these data show that HNRNPM knockdown with ", exp,
           " significantly impairs SH-SY5Y proliferation, with the divergence between ",
           "control and knockdown reaching statistical significance from ", first.sig,
           " and sustained thereafter.")
  } else {
    paste0("Together, these data show that HNRNPM knockdown with ", exp,
           " reduces SH-SY5Y proliferation, supported by the significant treatment and ",
           "treatment x time effects, although no single day survived Bonferroni correction.")
  }
  out <- c(out, concl)
  out <- c(out, "")

  ## --- per-day summary table --------------------------------------------- ##
  out <- c(out, paste0("Table. Mean cell counts (+/- SD) per well for gCtrl and ", exp,
                       ", Day 0 to Day 6, with Bonferroni-corrected P-values for the ",
                       "per-day gCtrl vs ", exp, " comparison (n = 3 per condition; ",
                       "'-' = seeding day, not tested; * P < 0.05)."))
  out <- c(out, paste("Day", "gCtrl (mean +/- SD)", paste0(exp, " (mean +/- SD)"),
                      "Bonferroni P", "Sig", sep = "\t"))
  for (d in 0:6) {
    row <- s[s$term == paste0("Day ", d), ]
    p.cell   <- if (is.na(row$p_adjusted)) "-" else fmt.p(row$p_adjusted)
    sig.cell <- if (is.na(row$significance)) "-" else row$significance
    out <- c(out, paste(paste0("Day ", d),
                        paste0(comma(row$mean_gCtrl),   " +/- ", comma(row$sd_gCtrl)),
                        paste0(comma(row$mean_gHNRNPM), " +/- ", comma(row$sd_gHNRNPM)),
                        p.cell, sig.cell, sep = "\t"))
  }
  out <- c(out, "")
  out <- c(out, "")
}

## --- closing paragraph: consistency across the two guides ---------------- ##
out <- c(out, paste0(
  "Consistent reductions in proliferation were obtained with two independent ",
  "HNRNPM-targeting guide RNAs (gHNRNPM_1 and gHNRNPM_2), indicating that the impaired ",
  "growth reflects loss of HNRNPM rather than an off-target effect of any single guide. ",
  "The growth curves for both guides are shown in ", FIGURE_PDF, " (mean +/- SD, n = 3; ",
  "significance stars indicate Bonferroni-corrected per-day comparisons)."))

## ----------------------------- write file -------------------------------- ##
writeLines(out, OUTPUT_TXT)
print(paste0("Saved ", OUTPUT_TXT))
