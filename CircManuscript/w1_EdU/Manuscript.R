# Manuscript.R
#
# Author: Juli Wang
#
# Purpose
# -------
# Generate a BMC Genomics-style Methods and Results text for the SH-SY5Y
# (dCas9-KRAB) HNRNPM-knockdown EdU S-phase assay in this folder. The Methods
# section reproduces the approved EdU methods draft, tailored to this dataset
# (two independent guide RNAs, and R used for the downstream statistics). The
# Results section is written live from the two stats files produced by the EdU
# analysis scripts (stats_studentttest.csv and stats_Welchttest.csv), so the
# reported means, SDs and p-values always match the current analysis, and it
# refers to the EdU bar-plot figures.
#
# Output: a plain-text file (OUTPUT_TXT) with a Methods section, a Results
# section reporting both t-tests, and a per-comparison summary table.
#
# Run this in RStudio (Source), or with: Rscript Manuscript.R

## ----------------------------- settings ---------------------------------- ##
STATS_STUDENT <- "stats_studentttest.csv"          # Student's t-test results
STATS_WELCH   <- "stats_Welchttest.csv"            # Welch's t-test results
FIG_STUDENT   <- "EdU_statistics_studentttest.pdf" # Student's bar-plot figure
FIG_WELCH     <- "EdU_statistics_Welchttest.pdf"   # Welch's bar-plot figure
OUTPUT_TXT    <- "EdU_Methods_Results.txt"          # the text file to write
RESULTS_DIR   <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w1_EdU/"
setwd(RESULTS_DIR)

guides <- c("gHNRNPM_1", "gHNRNPM_2")   # the two guide-RNA lines
dys    <- c("Day 1", "Day 4")           # the two harvest time points

## ----------------------------- helpers ----------------------------------- ##
# Format a p-value for prose: exact to 3 significant figures.
fmt.p <- function(p) paste0(signif(p, 3))
# Format a mean +/- SD percentage to one decimal place.
pct <- function(m, s) paste0(format(round(m, 1), nsmall = 1), " +/- ",
                             format(round(s, 1), nsmall = 1), "%")

## ----------------------------- read stats -------------------------------- ##
st <- read.csv(STATS_STUDENT, stringsAsFactors = FALSE, check.names = FALSE)
we <- read.csv(STATS_WELCH,   stringsAsFactors = FALSE, check.names = FALSE)

## ----------------------------- Methods text ------------------------------ ##
# Reproduced from the approved EdU methods draft
# (EdU_methods_draft_Carmen_JW_tracked_accepted.docx). Two edits were made to
# match this dataset: (1) two independent guide RNAs rather than one, and
# (2) the statistical-analysis paragraph now states that gating outputs were
# analysed in R (both Student's and Welch's t-tests), replacing the original
# GraphPad Prism / single-test wording.
out <- c()
out <- c(out, "Methods")
out <- c(out, "")
out <- c(out, "EdU proliferation and cell-cycle analysis")
out <- c(out, paste0(
  "Proliferation was assessed by 5-ethynyl-2′-deoxyuridine (EdU) incorporation using the ",
  "Click-iT Plus EdU Pacific Blue Flow Cytometry Assay Kit (Thermo Fisher Scientific), following ",
  "the manufacturer’s instructions with minor adaptations. SH-SY5Y cells stably expressing a ",
  "catalytically inactive Cas9 fused to a Krüppel-associated box repressor domain (dCas9-KRAB) ",
  "were used to generate control and knockdown lines. Two independent guide RNAs (gRNAs) targeting ",
  "HNRNPM (gHNRNPM_1 and gHNRNPM_2) were used to achieve CRISPR interference (CRISPRi)-mediated ",
  "knockdown, and each was compared with a line expressing a non-targeting control gRNA (gCtrl). ",
  "Cells were maintained in Dulbecco’s Modified Eagle Medium (DMEM) supplemented with 10% fetal ",
  "bovine serum (FBS) and seeded at equal density one or four days before the assay. For labelling, ",
  "EdU was added directly to the existing culture medium to a final concentration of 10 µM, and ",
  "cells were incubated for 2 hours at 37 °C. Matched no-EdU control wells received an equivalent ",
  "volume of EdU-free medium."))
out <- c(out, "")
out <- c(out, paste0(
  "Cells were washed with phosphate-buffered saline (PBS), dissociated with trypsin (~5 min, 37 °C), ",
  "neutralised in medium, and pelleted (250 g, 5 min). Pellets were washed in 1% bovine serum albumin ",
  "(BSA) in PBS, fixed in Click-iT fixative for 15 min at room temperature in the dark, washed again in ",
  "1% BSA in PBS, and permeabilised in 1× Click-iT saponin-based permeabilisation and wash reagent ",
  "for 15 min at room temperature. EdU was detected by adding the freshly prepared Click-iT Plus reaction ",
  "cocktail (containing the Pacific Blue picolyl azide) and incubating for 30 min at room temperature, ",
  "protected from light, followed by a wash in permeabilisation and wash reagent."))
out <- c(out, "")
out <- c(out, paste0(
  "DNA was counterstained by resuspending cells in FxCycle PI/RNase staining solution (Thermo Fisher ",
  "Scientific), which contains propidium iodide (PI) and ribonuclease (RNase), and incubating for 30 min ",
  "at room temperature in the dark. Immediately before acquisition, samples were passed through a ",
  "cell-strainer cap to remove aggregates and kept on ice."))
out <- c(out, "")
out <- c(out, paste0(
  "Samples were acquired on a BD LSRFortessa X20 flow cytometer running BD FACSDiva software v8.0.1, with ",
  "EdU–Pacific Blue detected in the violet-excited V450 channel and PI in the yellow-green-excited ",
  "YG-610 channel, and 10,000 events were recorded per sample. Detector voltages and analysis gates were ",
  "established using unstained cells, single-stained EdU-only cells, and single-stained PI-only cells. ",
  "Gating was applied hierarchically. Intact cells were first selected on forward-scatter area versus ",
  "side-scatter area (FSC-A versus SSC-A). Single cells were then resolved on forward-scatter height ",
  "versus forward-scatter area (FSC-H versus FSC-A), and cell-cycle populations were defined on a ",
  "bivariate plot of EdU–Pacific Blue (V450-A) versus DNA content (PI, YG-610-A), with polygonal gates ",
  "delineating G0/G1 (2N DNA, EdU-low), S phase (EdU-high), and G2/M (4N DNA, EdU-low). Flow cytometry ",
  "data were gated and quantified in BD FACSDiva v8.0.1."))
out <- c(out, "")
out <- c(out, "Statistical analysis")
out <- c(out, paste0(
  "The per-sample percentage of single cells in S phase was determined by hierarchical gating of the ",
  "flow-cytometry standard (.fcs) acquisition files in BD FACSDiva v8.0.1; these S-phase percentages ",
  "(not the raw .fcs files) were the input for statistical testing. Downstream statistical analysis and ",
  "data visualisation were performed in R (v4.3.3). For each guide RNA, the S-phase percentage was ",
  "compared between the gCtrl and gHNRNPM lines at each time point (Day 1 and Day 4) using an unpaired ",
  "two-tailed t-test, with n = 3 independent replicates per line; both the Student’s (equal-variance) ",
  "and Welch’s (unequal-variance) formulations are reported. Data are presented as mean ± standard ",
  "deviation (SD), and a difference was considered statistically significant at p < 0.05."))
out <- c(out, "")
out <- c(out, "")

## ----------------------------- Results text ------------------------------ ##
out <- c(out, "Results")
out <- c(out, "")
out <- c(out, "HNRNPM knockdown reduces S-phase entry in SH-SY5Y cells")
out <- c(out, "")
out <- c(out, paste0(
  "EdU incorporation was used to quantify the fraction of cells in S phase in control (gCtrl) and ",
  "HNRNPM-knockdown SH-SY5Y cells, using two independent guide RNAs (gHNRNPM_1 and gHNRNPM_2), one and ",
  "four days after plating (n = 3 per line; Fig. 1). For both guides and at both time points, HNRNPM ",
  "knockdown lowered the S-phase fraction relative to control."))
out <- c(out, "")

for (g in guides) {
  d1s <- st[st$guide == g & st$day == "Day 1", ]
  d4s <- st[st$guide == g & st$day == "Day 4", ]
  d1w <- we[we$guide == g & we$day == "Day 1", ]
  d4w <- we[we$guide == g & we$day == "Day 4", ]

  out <- c(out, paste0(
    "With ", g, ", the S-phase fraction was reduced from ", pct(d1s$mean_gCtrl, d1s$sd_gCtrl),
    " in gCtrl to ", pct(d1s$mean_gHNRNPM, d1s$sd_gHNRNPM), " at Day 1 (Student’s t-test P = ",
    fmt.p(d1s$p_value), ", ", d1s$significance, "; Welch’s t-test P = ", fmt.p(d1w$p_value),
    ", ", d1w$significance, "), and from ", pct(d4s$mean_gCtrl, d4s$sd_gCtrl), " to ",
    pct(d4s$mean_gHNRNPM, d4s$sd_gHNRNPM), " at Day 4 (Student’s P = ", fmt.p(d4s$p_value),
    ", ", d4s$significance, "; Welch’s P = ", fmt.p(d4w$p_value), ", ", d4w$significance, ")."))
  out <- c(out, "")
}

out <- c(out, paste0(
  "Together, EdU labelling shows that HNRNPM knockdown reduces S-phase entry in SH-SY5Y cells with two ",
  "independent guide RNAs, consistent with the slower proliferation seen in the growth-curve assay. The ",
  "reduction was statistically significant at both time points by both tests for gHNRNPM_2; for ",
  "gHNRNPM_1 it was significant at both days by the Student’s t-test and at Day 4 by the more ",
  "conservative Welch’s t-test (Day 1 being borderline, reflecting the very small control-group ",
  "variance and n = 3). The corresponding bar plots are provided in ", FIG_STUDENT, " (Student’s) and ",
  FIG_WELCH, " (Welch’s), showing mean ± SD with individual replicates and the per-comparison p-value."))
out <- c(out, "")
out <- c(out, "")

## ----------------------------- summary table ----------------------------- ##
out <- c(out, paste0("Table. Percentage of single cells in S phase (mean +/- SD) for gCtrl and each ",
                     "HNRNPM-knockdown line at Day 1 and Day 4, with unpaired two-tailed t-test p-values ",
                     "(n = 3 per line; * P < 0.05, ** P < 0.01, *** P < 0.001; ns = not significant)."))
out <- c(out, paste("Guide", "Day", "gCtrl (%S)", "gHNRNPM (%S)", "n",
                    "Student P (sig)", "Welch P (sig)", sep = "\t"))
for (g in guides) {
  for (dy in dys) {
    rs <- st[st$guide == g & st$day == dy, ]
    rw <- we[we$guide == g & we$day == dy, ]
    out <- c(out, paste(g, dy,
                        pct(rs$mean_gCtrl, rs$sd_gCtrl),
                        pct(rs$mean_gHNRNPM, rs$sd_gHNRNPM),
                        rs$n_ctrl,
                        paste0(fmt.p(rs$p_value), " (", rs$significance, ")"),
                        paste0(fmt.p(rw$p_value), " (", rw$significance, ")"),
                        sep = "\t"))
  }
}

## ----------------------------- write file -------------------------------- ##
writeLines(out, OUTPUT_TXT, useBytes = TRUE)
print(paste0("Saved ", OUTPUT_TXT))
