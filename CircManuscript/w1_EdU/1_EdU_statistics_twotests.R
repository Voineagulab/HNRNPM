#!/usr/bin/env Rscript
###############################################################################
# EdU_statistics.R  (two-tests variant)
#
# Author: Juli Wang
#
# Purpose
# -------
# Analyse the EdU S-phase data in `Raw_percentage.xlsx` (sheet "% S-phase").
# SH-SY5Y cells stably expressing dCas9-KRAB were transduced with either a
# non-targeting control gRNA (gCtrl) or one of two HNRNPM-targeting gRNAs
# (gHNRNPM_1, gHNRNPM_2), selected with puromycin, then plated and harvested
# either 1 day (Day 1) or 4 days (Day 4) post-plating. The measurement is the
# percentage of cells in the S-phase gate.
#
# For each guide (gHNRNPM_1, gHNRNPM_2) and each time point (Day 1, Day 4) the
# script tests whether the % of cells in S-phase differs between the
# HNRNPM-targeting guide and the matched gCtrl (4 comparisons total), then
# draws an annotated bar plot for every comparison. Each comparison is a
# single-column-width figure (85 mm, BMC format), and all four are written
# one-per-page into a single multi-page PDF (4 pages).
#
# Statistical tests (this variant reports BOTH)
# ---------------------------------------------
# Two-sample, two-sided t-tests comparing gCtrl vs gHNRNPM for each guide/day.
#   * The p-value printed in the bracket ABOVE the bars is the STUDENT's
#     (equal-variance) t-test p-value  -> identical to the value shown in
#     EdU_statistics_studentttest.pdf.
#   * A CAPTION beneath each panel additionally reports the WELCH's
#     (unequal-variance) t-test p-value.
# Both p-values use the same format: an italicised p in scientific notation
# with two decimal places, e.g. p = 5.46 x 10^-3.
#
# Usage
# -----
#   Rscript 1_EdU_statistics_twotests.R
# (run from the folder containing Raw_percentage.xlsx, or edit the paths below)
###############################################################################

## ----------------------------- settings ---------------------------------- ##
INPUT_XLSX  <- "Raw_percentage.xlsx"      # input workbook
SHEET_NAME  <- "% S-phase"                # sheet with the raw percentages
ALTERNATIVE <- "two.sided"                # two-sided: test for any change
OUTPUT_PDF  <- "EdU_statistics_twotests.pdf"   # output figure
OUTPUT_CSV  <- "stats_twotests.csv"            # both t-tests' results
RESULTS_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w1_EdU/"
setwd(RESULTS_DIR)

## Figure size (BMC single-column width). Units are millimetres; converted to
## inches for the pdf() device below (1 inch = 25.4 mm).
FIG_WIDTH_MM  <- 85    # BMC single-column width
FIG_HEIGHT_MM <- 95    # keep well under the 225 mm max height
BAR_WIDTH     <- 0.45  # narrower bars (0-1); lower = thinner bars, more gap
BAR_SPACING   <- 0.6   # centre-to-centre distance between the two bars
                       # (smaller = bars closer together; must be > BAR_WIDTH)
BASE_FONT     <- 12    # minimum font size (pt) for all text in the plots

## Resolve paths relative to this script's directory when possible, so the
## script works whether it is sourced or run with Rscript from another folder.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(normalizePath(dirname(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(dirname(sys.frames()[[1]]$ofile)))
  getwd()
}
script_dir <- tryCatch(get_script_dir(), error = function(e) getwd())
resolve <- function(p) if (file.exists(p)) p else file.path(script_dir, p)
INPUT_XLSX <- resolve(INPUT_XLSX)
if (!file.exists(INPUT_XLSX))
  stop("Cannot find input workbook: ", INPUT_XLSX)
OUTPUT_PDF <- if (dirname(OUTPUT_PDF) == ".") file.path(script_dir, OUTPUT_PDF) else OUTPUT_PDF

## ------------------------- package management ----------------------------- ##
need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
invisible(lapply(c("readxl", "ggplot2"), need))
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

## --------------------------- read the data -------------------------------- ##
# The sheet has two stacked blocks with the following fixed layout
# (1-based sheet rows/columns; empty spacer column C separates Day1 and Day4):
#
#   Block gHNRNPM_1   Day1: col A = gCtrl, col B = gHNRNPM_1   (rows 5-7)
#                     Day4: col D = gCtrl, col E = gHNRNPM_1   (rows 5-7)
#   Block gHNRNPM_2   Day1: col A = gCtrl, col B = gHNRNPM_2   (rows 12-14)
#                     Day4: col D = gCtrl, col E = gHNRNPM_2   (rows 12-14)
raw <- as.data.frame(read_excel(INPUT_XLSX, sheet = SHEET_NAME,
                                col_names = FALSE, .name_repair = "minimal"))

getvals <- function(rows, col) {
  v <- suppressWarnings(as.numeric(raw[rows, col]))
  v[!is.na(v)]
}

comparisons <- list(
  list(guide = "gHNRNPM_1", day = "Day 1",
       ctrl = getvals(5:7, 1),   treat = getvals(5:7, 2)),
  list(guide = "gHNRNPM_1", day = "Day 4",
       ctrl = getvals(5:7, 4),   treat = getvals(5:7, 5)),
  list(guide = "gHNRNPM_2", day = "Day 1",
       ctrl = getvals(12:14, 1), treat = getvals(12:14, 2)),
  list(guide = "gHNRNPM_2", day = "Day 4",
       ctrl = getvals(12:14, 4), treat = getvals(12:14, 5))
)

## ------------------------------ helpers ----------------------------------- ##
sem   <- function(x) sd(x) / sqrt(length(x))
stars <- function(p) {
  if (is.na(p))      "n/a"
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 5e-2) "*"
  else               "ns"
}
# Format a p-value as scientific notation with two decimal places, returned as
# a plotmath string (rendered with parse = TRUE), e.g. 0.00546 -> 5.46 x 10^-3.
# The mantissa is quoted so the two decimals (incl. trailing zeros) are kept.
p_expr <- function(p) {
  if (is.na(p)) return("italic(p)*' = NA'")
  exponent <- floor(log10(p))
  mantissa <- round(p / 10^exponent, 2)
  if (mantissa >= 10) { mantissa <- mantissa / 10; exponent <- exponent + 1 }
  paste0("italic(p)=='", formatC(mantissa, format = "f", digits = 2),
         "'%*%10^", exponent)
}

## ------------------- run tests + build the four panels -------------------- ##
plots        <- list()
summary_rows <- list()

for (i in seq_along(comparisons)) {
  cmp   <- comparisons[[i]]
  ctrl  <- cmp$ctrl
  treat <- cmp$treat

  # Both tests on the same data: Student's (equal var) and Welch's (unequal var).
  tt_s <- t.test(ctrl, treat, alternative = ALTERNATIVE, var.equal = TRUE)
  tt_w <- t.test(ctrl, treat, alternative = ALTERNATIVE, var.equal = FALSE)
  p_s  <- tt_s$p.value    # Student's -> shown in the bracket above the bars
  p_w  <- tt_w$p.value    # Welch's   -> shown in the caption

  # tidy frames for plotting. Bars are placed at explicit numeric x-positions
  # (xpos) that are BAR_SPACING apart, so the two bars can be drawn close
  # together regardless of the panel width.
  grp_levels <- c("gCtrl", cmp$guide)
  xpos       <- c(1, 1 + BAR_SPACING)          # gCtrl, gHNRNPM
  raw_df <- data.frame(
    group = factor(rep(grp_levels, c(length(ctrl), length(treat))),
                   levels = grp_levels),
    x     = rep(xpos, c(length(ctrl), length(treat))),
    value = c(ctrl, treat)
  )
  summ_df <- data.frame(
    group = factor(grp_levels, levels = grp_levels),
    x     = xpos,
    mean  = c(mean(ctrl), mean(treat)),
    sem   = c(sem(ctrl),  sem(treat))
  )

  # y-position for the significance bracket
  ymax   <- max(raw_df$value, summ_df$mean + summ_df$sem)
  bar_y  <- ymax * 1.10
  tick   <- ymax * 0.02
  text_y <- bar_y + ymax * 0.06

  cols <- c("#9E9E9E", "#D55E00"); names(cols) <- grp_levels
  ann_size <- BASE_FONT / .pt                  # convert pt -> ggplot size units

  # Welch caption (same italic-p scientific format as the bracket p).
  welch_cap <- parse(text = paste0("\"Welch's two-sample t-test: \"*", p_expr(p_w)))

  pg <- ggplot(summ_df, aes(x = x, y = mean, fill = group)) +
    geom_col(width = BAR_WIDTH, colour = "black") +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = BAR_WIDTH * 0.4) +
    geom_jitter(data = raw_df, aes(x = x, y = value),
                width = BAR_WIDTH * 0.25, height = 0, size = 2,
                shape = 21, fill = "white", colour = "black",
                inherit.aes = FALSE) +
    # significance bracket (Student's t-test p-value)
    annotate("segment", x = xpos[1], xend = xpos[2], y = bar_y, yend = bar_y) +
    annotate("segment", x = xpos[1], xend = xpos[1], y = bar_y, yend = bar_y - tick) +
    annotate("segment", x = xpos[2], xend = xpos[2], y = bar_y, yend = bar_y - tick) +
    annotate("text", x = mean(xpos), y = text_y,
             label = paste0("'", stars(p_s), "   '*", p_expr(p_s)), parse = TRUE,
             size = ann_size, fontface = "bold") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(breaks = xpos, labels = grp_levels,
                       expand = expansion(add = 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = paste0(cmp$guide, " vs gCtrl - ", cmp$day),
         subtitle = "Student's two-sample t-test",
         x = NULL, y = "% cells in S-phase",
         caption = welch_cap) +
    # BASE_FONT (12 pt) is the floor; title/subtitle sit at or above it
    theme_classic(base_size = BASE_FONT) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = BASE_FONT),
          plot.subtitle = element_text(hjust = 0.5, size = BASE_FONT),
          plot.caption  = element_text(hjust = 0.5, size = BASE_FONT - 2, colour = "black"),
          axis.text     = element_text(size = BASE_FONT),
          axis.title    = element_text(size = BASE_FONT),
          axis.text.x   = element_text(face = "bold"))

  plots[[i]] <- pg

  summary_rows[[i]] <- data.frame(
    guide        = cmp$guide,
    day          = cmp$day,
    n_ctrl       = length(ctrl),
    n_treat      = length(treat),
    mean_gCtrl   = round(mean(ctrl), 2),
    sd_gCtrl     = round(sd(ctrl),   2),
    mean_gHNRNPM = round(mean(treat),2),
    sd_gHNRNPM   = round(sd(treat),  2),
    student_t    = round(unname(tt_s$statistic), 3),
    student_df   = round(unname(tt_s$parameter), 2),
    student_p    = signif(p_s, 4),
    student_sig  = stars(p_s),
    welch_t      = round(unname(tt_w$statistic), 3),
    welch_df     = round(unname(tt_w$parameter), 2),
    welch_p      = signif(p_w, 4),
    welch_sig    = stars(p_w),
    stringsAsFactors = FALSE
  )
}

## ---------------------- print results to console -------------------------- ##
results <- do.call(rbind, summary_rows)
cat("\n=== EdU % S-phase: gHNRNPM vs gCtrl (Student's + Welch's) ===\n\n")
print(results, row.names = FALSE)
cat("\n")

## --------------------- write all test results to CSV ---------------------- ##
# Both t-tests for the four gCtrl vs gHNRNPM comparisons, with group means,
# SDs, n, and each test's t-statistic, df, p-value and significance.
write.csv(results, OUTPUT_CSV, row.names = FALSE)
cat("Saved test results to: ", OUTPUT_CSV, "\n", sep = "")

## ------------------------ write the PDF figure ---------------------------- ##
# One comparison per page at BMC single-column width -> a 4-page PDF.
# onefile = TRUE (default) keeps every page in the same file.
pdf(OUTPUT_PDF,
    width  = FIG_WIDTH_MM  / 25.4,
    height = FIG_HEIGHT_MM / 25.4,
    onefile = TRUE)
for (pg in plots) print(pg)
invisible(dev.off())

cat("Saved ", length(plots), "-page figure to: ", OUTPUT_PDF, "\n", sep = "")
