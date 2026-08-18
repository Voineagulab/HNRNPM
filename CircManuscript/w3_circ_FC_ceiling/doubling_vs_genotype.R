# doubling_vs_genotype.R
#
# Author: Juli Wang
#
# Purpose
# -------
# Ask whether the mild proliferation slow-down caused by HNRNPM knockdown could,
# by itself, account for the increased circRNA levels in the knockdown. The model
# is the steady-state single-compartment dilution of a stable RNA: for the same
# transcription, a slower-dividing cell dilutes a stable RNA less, so its per-cell
# level rises. The proliferation-only fold change (KD / control) is
#
#     FC(half_life) = (k_deg + mu_ctrl) / (k_deg + mu_KD)
#
# with k_deg = ln2 / half_life, mu = ln2 / doubling_time, and
# doubling_time_KD = R * doubling_time_ctrl. R is the KD-to-control doubling-time
# ratio. As half_life -> infinity the fold change approaches its CEILING = R
# (the largest circRNA fold change that reduced division could ever produce).
#
# R is estimated two ways per guide:
#   * EdU-derived  : R = %S_ctrl / %S_KD  (assumes constant S-phase duration),
#                    reported at Day 1 and Day 4.
#   * Growth-curve : R = Td_KD / Td_ctrl, from a log-linear fit of the cell
#                    counts restricted to the exponential phase. The exponential
#                    window is chosen OBJECTIVELY, not by eye: among all
#                    consecutive-day windows of >= MIN_WIN timepoints, the one
#                    whose log-linear fit is most linear (highest mean R^2 across
#                    the two treatments) is kept, and the SAME window is used for
#                    gCtrl and gHNRNPM so that R is a fair ratio. This avoids the
#                    plateau (D5-D6, as cultures approach confluence) biasing the
#                    doubling times. The full per-window scan is written to
#                    growth_window_selection.csv for transparency.
#
# The plotted ceilings are compared with the SMALLEST fold change among the 36
# HNRNPM-independent circRNAs (Supplementary Table S5B, column logFC_BSJ): if even
# that smallest circRNA fold change exceeds the ceiling, reduced proliferation
# cannot fully explain the circRNA increases.
#
# Outputs:
#   circ_FC_ceiling.pdf       - 2 pages: Main (Day 1), Supplementary (Day 4)
#   circ_FC_ceiling.pptx      - 2 slides: Main (Day 1), Supplementary (Day 4)
#   growth_window_selection.csv - per-window R^2 scan + the chosen window
#   manuscript_inputs.csv     - every quantity Manuscript.R needs, so the text
#                               and the figure share ONE source of truth. Run
#                               this script before Manuscript.R.
#
# Run in RStudio (Source), or: Rscript doubling_vs_genotype.R

## ----------------------------- settings ---------------------------------- ##
RESULTS_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w3_circ_FC_ceiling/"
EDU_STATS   <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w1_EdU/stats_studentttest.csv"
GROWTH_XLSX <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w2_growth_curve/growth_curve.xlsx"
CIRC_XLSX   <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/Supplementary_Table/Supplementary_Tables_JW_2.xlsx"
CIRC_SHEET  <- "S5B_independent_circRNAs"   # header is on row 2 (row 1 is a title) -> skip = 1

OUTPUT_PDF    <- "circ_FC_ceiling.pdf"
OUTPUT_PPTX   <- "circ_FC_ceiling.pptx"
OUTPUT_WIN    <- "growth_window_selection.csv"   # per-window R^2 scan + chosen window
OUTPUT_INPUTS <- "manuscript_inputs.csv"         # every quantity Manuscript.R needs

setwd(RESULTS_DIR)

guides   <- c("gHNRNPM_1", "gHNRNPM_2")
gc_sheet <- c(gHNRNPM_1 = "GrowthCurve_gHNRNPM_1_raw",
              gHNRNPM_2 = "GrowthCurve_gHNRNPM_2_raw")
days     <- 0:6

## ------------------- Exponential-window selection settings ---------------- ##
## MIN_WIN = the minimum number of TIMEPOINTS (distinct days) a candidate window
##           must contain to be considered for the log-linear doubling-time fit.
##
## Why a floor at all?
##   Cell counts flatten toward D5-D6 as cultures approach confluence, so a fit
##   over the whole D0-D6 range underestimates the true exponential growth rate.
##   We therefore search for the exponential window objectively (make_windows +
##   the R^2 scan below) instead of choosing it by eye. But raw R^2 tends to rise
##   mechanically as a window shrinks (fewer points left to deviate from the
##   line), so with no floor the search would drift toward trivially short
##   windows that "win" on R^2 without representing a real growth phase.
##
## Why MIN_WIN = 4 specifically (not 3):
##   1. Curvature needs >=4 points to be visible. One stated selection criterion
##      is "shrink until the residuals stop showing curvature"; with only 3
##      timepoints a straight line always fits with a single residual and there
##      is no curvature signal to test. 4 is the smallest window that can.
##   2. It limits the short-window R^2 bias described above.
##   3. Per-replicate doubling times (Td_reps) use only the window's timepoints;
##      4 points per replicate keeps those slope estimates reasonably stable
##      (the KD lines are already the noisier fits, R^2 ~0.77-0.83).
##
## Empirical check (run on this dataset, both guides):
##   Lowering MIN_WIN to 3 changes NOTHING here. The newly-allowed 3-timepoint
##   windows (D0-D2, D1-D3, ...) never beat D0-D3 / D0-D4 on mean R^2, so the
##   selected windows -- and therefore the doubling times -- are identical:
##       gHNRNPM_1: window D0-D3, Td_ctrl = 48.5 h, Td_KD = 73.5 h, R = 1.513
##       gHNRNPM_2: window D0-D4, Td_ctrl = 59.6 h, Td_KD = 90.2 h, R = 1.514
##   i.e. D0-D2 is not more log-linear than D0-D3 in this data. We keep MIN_WIN =
##   4 as the more conservative, curvature-aware default; it costs nothing here.
##
## USE_ADJ_R2 = FALSE ranks windows by mean raw R^2 across the two treatments.
##   Set TRUE to rank by mean ADJUSTED R^2 instead (penalises extra points more,
##   an alternative guard against the short-window bias). Ties in the ranking are
##   broken toward the LONGER window, then the earlier start (see growth_fit).
MIN_WIN     <- 4
USE_ADJ_R2  <- FALSE

## ----------------------------- packages ---------------------------------- ##
library(readxl)     # read the .xlsx inputs
library(ggplot2)    # plotting
library(patchwork)  # combine the two guide panels
library(officer)    # write the .pptx (figure embedded as a high-resolution image)

ln2 <- log(2)

## ----------------------------- inputs ------------------------------------ ##
## EdU S-phase means -> EdU-derived doubling-time ratio R = %S_ctrl / %S_KD
edu <- read.csv(EDU_STATS, stringsAsFactors = FALSE, check.names = FALSE)
R_edu <- function(guide, day) {
  r <- edu[edu$guide == guide & edu$day == day, ]
  r$mean_gCtrl / r$mean_gHNRNPM
}

## Growth-curve doubling times on the objectively-chosen exponential window.
## R = Td_KD / Td_ctrl, plus control Td (hours). See header + settings above.

## pooled log-linear fit of log(count) ~ day for one treatment over a day-set;
## returns slope (per day), R^2 and adjusted R^2 (reps treated as points)
fit_pool <- function(rows, day_set) {
  cnt <- c(); day <- c()
  for (r in 1:nrow(rows)) {
    for (d in day_set) { cnt <- c(cnt, rows[r, paste0("D", d)]); day <- c(day, d) }
  }
  m  <- lm(log(cnt) ~ day)
  sm <- summary(m)
  c(slope = unname(coef(m)[2]), r2 = sm$r.squared, adj_r2 = sm$adj.r.squared)
}

## per-replicate doubling times (hours) over a day-set -> mean/SD for uncertainty
Td_reps <- function(rows, day_set) {
  td <- c()
  for (r in 1:nrow(rows)) {
    cnt <- as.numeric(rows[r, paste0("D", day_set)])
    td  <- c(td, (ln2 / unname(coef(lm(log(cnt) ~ day_set))[2])) * 24)
  }
  td
}

## all consecutive-day windows with >= min_win timepoints
make_windows <- function(days, min_win) {
  w <- list()
  for (s in seq_along(days)) for (e in seq_along(days)) {
    if (e >= s && (e - s + 1) >= min_win) w[[length(w) + 1]] <- days[s:e]
  }
  w
}

growth_fit <- function(sheet) {
  wide <- as.data.frame(read_excel(GROWTH_XLSX, sheet = sheet))
  rc <- wide[wide$Cell_line == "gCtrl",   ]
  rk <- wide[wide$Cell_line == "gHNRNPM", ]

  ## score every candidate window by mean linearity across the two treatments
  score_col <- if (USE_ADJ_R2) "adj_r2" else "r2"
  scan <- data.frame()
  for (w in make_windows(days, MIN_WIN)) {
    fc <- fit_pool(rc, w); fk <- fit_pool(rk, w)
    scan <- rbind(scan, data.frame(
      start    = min(w), end = max(w), n_days = length(w),
      r2_ctrl  = unname(fc["r2"]),     r2_kd  = unname(fk["r2"]),
      adjr2_ctrl = unname(fc["adj_r2"]), adjr2_kd = unname(fk["adj_r2"]),
      Td_ctrl_h = (ln2 / unname(fc["slope"])) * 24,
      Td_kd_h   = (ln2 / unname(fk["slope"])) * 24,
      score     = (unname(fc[score_col]) + unname(fk[score_col])) / 2))
  }
  scan$R_growth <- scan$Td_kd_h / scan$Td_ctrl_h   # R = Td_KD / Td_ctrl

  ## best window: highest mean (adj)R^2; ties -> longer window, then earlier start
  scan <- scan[order(-scan$score, -scan$n_days, scan$start), ]
  scan$chosen <- FALSE; scan$chosen[1] <- TRUE
  best <- scan[1, ]
  win  <- best$start:best$end

  ## per-replicate doubling times on the chosen window (for mean +/- SD)
  tdc <- Td_reps(rc, win); tdk <- Td_reps(rk, win)

  list(
    win        = win,
    win_label  = paste0("D", best$start, "-D", best$end),
    scan       = scan,
    Td_ctrl_h  = best$Td_ctrl_h,
    Td_kd_h    = best$Td_kd_h,
    R_growth   = best$R_growth,
    r2_ctrl    = best$r2_ctrl,
    r2_kd      = best$r2_kd,
    Td_ctrl_h_repmean = mean(tdc), Td_ctrl_h_repsd = sd(tdc),
    Td_kd_h_repmean   = mean(tdk), Td_kd_h_repsd   = sd(tdk),
    R_growth_repmean  = mean(tdk) / mean(tdc))
}
gc <- lapply(gc_sheet, growth_fit)           # per guide

## write the full per-window scan (both guides) for transparency
scan_all <- do.call(rbind, lapply(names(gc), function(g) {
  s <- gc[[g]]$scan; s$guide <- g
  s[, c("guide", setdiff(names(s), "guide"))]
}))
write.csv(scan_all, OUTPUT_WIN, row.names = FALSE)
print(paste0("Saved ", OUTPUT_WIN))

## Fold changes among the independent circRNAs (S5B, logFC_BSJ)
circ <- as.data.frame(read_excel(CIRC_XLSX, sheet = CIRC_SHEET, skip = 1))
logfc <- circ$logFC_BSJ[!is.na(circ$logFC_BSJ)]
FC_min <- 2 ^ min(logfc)                     # smallest observed circRNA fold change
FC_max <- 2 ^ max(logfc)                     # largest observed circRNA fold change
n_circ <- length(logfc)

## -------- export every quantity Manuscript.R needs, as one tidy CSV ------- ##
## Manuscript.R reads THIS file rather than re-deriving anything, so the growth-
## fit, EdU and circRNA logic lives in exactly one place (this script). The file
## is tidy (metric, guide, value); guide = "all" for whole-dataset scalars.
## Run this script before Manuscript.R.
# human-readable explanation of each metric (filled into the "Note" column)
notes <- c(
  win_label  = "The time window with the highest exponential growth score",
  Td_ctrl_h  = "Doubling time of gCtrl (hours)",
  Td_kd_h    = "Doubling time of gHNRNPM knockdown (hours)",
  R_growth   = "Ratio of doubling time (KD/Ctrl)",
  r2_ctrl    = "R-squared of the gCtrl log-linear fit over the chosen window",
  r2_kd      = "R-squared of the gHNRNPM log-linear fit over the chosen window",
  R_edu_Day1 = "EdU-derived doubling-time ratio at Day 1 (%S_ctrl / %S_KD)",
  R_edu_Day4 = "EdU-derived doubling-time ratio at Day 4 (%S_ctrl / %S_KD)",
  FC_min     = "Smallest fold change among the independent circRNAs (logFC_BSJ)",
  FC_max     = "Largest fold change among the independent circRNAs (logFC_BSJ)",
  n_circ     = "Number of independent circRNAs")

minp <- data.frame(metric = character(), guide = character(),
                   value = character(), Note = character(), stringsAsFactors = FALSE)
addv <- function(metric, guide, value)
  minp <<- rbind(minp, data.frame(metric = metric, guide = guide,
                                  value = as.character(value),
                                  Note = unname(notes[metric]), stringsAsFactors = FALSE))
for (g in guides) {
  addv("win_label",  g, gc[[g]]$win_label)   # chosen exponential window, e.g. "D0-D3"
  addv("Td_ctrl_h",  g, gc[[g]]$Td_ctrl_h)   # control doubling time (hours)
  addv("Td_kd_h",    g, gc[[g]]$Td_kd_h)     # knockdown doubling time (hours)
  addv("R_growth",   g, gc[[g]]$R_growth)    # growth-curve doubling-time ratio Td_KD/Td_ctrl
  addv("r2_ctrl",    g, gc[[g]]$r2_ctrl)     # log-linear fit R^2, control, chosen window
  addv("r2_kd",      g, gc[[g]]$r2_kd)       # log-linear fit R^2, knockdown, chosen window
  addv("R_edu_Day1", g, R_edu(g, "Day 1"))   # EdU-derived ratio, Day 1
  addv("R_edu_Day4", g, R_edu(g, "Day 4"))   # EdU-derived ratio, Day 4
}
addv("FC_min", "all", FC_min)                # smallest independent-circRNA fold change
addv("FC_max", "all", FC_max)                # largest independent-circRNA fold change
addv("n_circ", "all", n_circ)                # number of independent circRNAs
write.csv(minp, OUTPUT_INPUTS, row.names = FALSE)
print(paste0("Saved ", OUTPUT_INPUTS))

## ----------------------------- model ------------------------------------- ##
prolif_fc <- function(half_life, T_ctrl_h, R) {
  k  <- ln2 / half_life
  mW <- ln2 / T_ctrl_h
  mK <- ln2 / (R * T_ctrl_h)
  (k + mW) / (k + mK)
}

# colours: EdU ceiling, growth-curve ceiling, circRNA reference line
col_edu    <- "#0072B2"
col_growth <- "#D55E00"
col_circ   <- "#111111"

## one panel: proliferation-only FC vs half-life, for one guide at one EdU day
panel <- function(guide, day) {
  Redu <- R_edu(guide, day)
  Rgro <- gc[[guide]]$R_growth
  Tctl <- gc[[guide]]$Td_ctrl_h
  win_lab <- gc[[guide]]$win_label

  half <- 10 ^ seq(log10(3), log10(1000), length.out = 400)
  df <- rbind(
    data.frame(half = half, fc = prolif_fc(half, Tctl, Redu), src = "EdU S-phase ratio"),
    data.frame(half = half, fc = prolif_fc(half, Tctl, Rgro), src = "Growth-curve doubling time"))

  y.top <- max(FC_min, Redu, Rgro) * 1.06

  ggplot(df, aes(half, fc, colour = src)) +
    # ceilings (dashed) and the smallest-circRNA-FC reference (solid black)
    geom_hline(yintercept = Redu, linetype = "dashed", colour = col_edu,    linewidth = 0.6) +
    geom_hline(yintercept = Rgro, linetype = "dashed", colour = col_growth, linewidth = 0.6) +
    geom_hline(yintercept = FC_min, colour = col_circ, linetype = "dotted", linewidth = 0.7) +
    geom_line(linewidth = 1) +
    annotate("text", x = 3, y = FC_min, hjust = 0, vjust = -0.5, size = 3, colour = col_circ,
             label = paste0("smallest independent circRNA FC = ", round(FC_min, 2),
                            "  (n = ", n_circ, ")")) +
    annotate("text", x = 1000, y = Redu, hjust = 1, vjust = -0.6, size = 3, colour = col_edu,
             label = paste0("EdU ceiling ratio = ", round(Redu, 2))) +
    annotate("text", x = 1000, y = Rgro, hjust = 1, vjust = -0.6, size = 3, colour = col_growth,
             label = paste0("growth-curve ceiling ratio = ", round(Rgro, 2))) +
    scale_colour_manual(values = c("EdU S-phase ratio" = col_edu,
                                   "Growth-curve doubling time" = col_growth),
                        name = NULL) +
    scale_x_log10(breaks = c(3, 10, 30, 100, 300, 1000),
                  labels = c("3", "10", "30", "100", "300", "1000")) +
    coord_cartesian(ylim = c(1.0, y.top)) +
    labs(title = paste0(guide, "  (", day, " S-phase)"),
         subtitle = paste0("growth-curve exponential window ", win_lab),
         x = "circRNA half-life (hours)",
         y = "circRNA fold increase (KD/Ctrl) attributable to\nreduced proliferation in HNRNPM-KD") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

## ----------------------------- figures ----------------------------------- ##
main.fig <- (panel("gHNRNPM_1", "Day 1") | panel("gHNRNPM_2", "Day 1")) +
  plot_annotation(title = "Main: proliferation-only circRNA fold-change ceiling (Day 1)")
supp.fig <- (panel("gHNRNPM_1", "Day 4") | panel("gHNRNPM_2", "Day 4")) +
  plot_annotation(title = "Supplementary: proliferation-only circRNA fold-change ceiling (Day 4)")

## ----------------------------- write PDF (2 pages) ----------------------- ##
pdf(OUTPUT_PDF, width = 12, height = 5.2)
print(main.fig)
print(supp.fig)
invisible(dev.off())
print(paste0("Saved ", OUTPUT_PDF))

## ----------------------------- write PPTX (2 slides) --------------------- ##
# Embed the same two figures as high-resolution images (no extra package needed).
png_main <- tempfile(fileext = ".png")
png_supp <- tempfile(fileext = ".png")
ggsave(png_main, main.fig, width = 12, height = 5.2, dpi = 200)
ggsave(png_supp, supp.fig, width = 12, height = 5.2, dpi = 200)

ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
ppt <- ph_with(ppt, external_img(png_main, width = 12, height = 5.2),
               location = ph_location(left = 0.15, top = 0.7, width = 9.7, height = 4.2))
ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
ppt <- ph_with(ppt, external_img(png_supp, width = 12, height = 5.2),
               location = ph_location(left = 0.15, top = 0.7, width = 9.7, height = 4.2))
print(ppt, target = OUTPUT_PPTX)
print(paste0("Saved ", OUTPUT_PPTX))

## ----------------------------- console summary --------------------------- ##
for (g in guides) {
  print(paste0(g,
    ": R_EdU(Day1) = ", round(R_edu(g, "Day 1"), 3),
    ", R_EdU(Day4) = ", round(R_edu(g, "Day 4"), 3),
    ", R_growth = ",   round(gc[[g]]$R_growth, 3),
    " [window ", gc[[g]]$win_label,
    ", R2 ctrl ", round(gc[[g]]$r2_ctrl, 3), " / KD ", round(gc[[g]]$r2_kd, 3), "]",
    ", Td_ctrl = ",    round(gc[[g]]$Td_ctrl_h, 1), " h",
    ", Td_KD = ",      round(gc[[g]]$Td_kd_h, 1), " h"))
  print(paste0("    per-replicate Td (chosen window): ctrl ",
    round(gc[[g]]$Td_ctrl_h_repmean, 1), " +/- ", round(gc[[g]]$Td_ctrl_h_repsd, 1),
    " h, KD ", round(gc[[g]]$Td_kd_h_repmean, 1), " +/- ", round(gc[[g]]$Td_kd_h_repsd, 1),
    " h, R (rep-mean) = ", round(gc[[g]]$R_growth_repmean, 3)))
}
print(paste0("Smallest independent circRNA FC = ", round(FC_min, 3), " (n = ", n_circ, ")"))
