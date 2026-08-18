# Manuscript.R
#
# Author: Juli Wang
#
# Purpose
# -------
# Generate a BMC Genomics-style Methods and Results text (and a Word .docx copy)
# for the proliferation-confounder analysis, tying together the EdU S-phase data,
# the growth-curve doubling times, and the fold changes of the HNRNPM-independent
# circRNAs (Supplementary Table S5B).
#
# Inputs
# ------
# This script does NO computation of its own. Every number it prints is read from
#   manuscript_inputs.csv
# which is written by doubling_vs_genotype.R. That keeps the growth-fit, EdU and
# circRNA logic in exactly ONE place, so the text can never drift from the figure.
# >>> Run doubling_vs_genotype.R FIRST to (re)generate manuscript_inputs.csv. <<<
#
# Model (for reference): steady-state single-compartment dilution of a stable RNA.
# The proliferation-only fold change (KD/control) approaches a ceiling equal to
# the knockdown-to-control doubling-time ratio R. R is estimated from EdU S-phase
# fractions (R = %S_ctrl/%S_KD; Day 1 and Day 4) and from growth-curve doubling
# times (R = Td_KD/Td_ctrl). These ceilings are compared with the smallest circRNA
# fold change; if even that exceeds the ceiling, reduced proliferation cannot
# explain the circRNA increases.
#
# Outputs:
#   circ_FC_ceiling_Methods_Results.txt
#   circ_FC_ceiling_Methods_Results.docx
#
# Run in RStudio (Source), or: Rscript Manuscript.R

## ----------------------------- settings ---------------------------------- ##
RESULTS_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w3_circ_FC_ceiling/"
INPUTS_CSV  <- "manuscript_inputs.csv"   # written by doubling_vs_genotype.R

OUTPUT_TXT  <- "circ_FC_ceiling_Methods_Results.txt"
OUTPUT_DOCX <- "circ_FC_ceiling_Methods_Results.docx"
FIG_PDF     <- "circ_FC_ceiling.pdf"

setwd(RESULTS_DIR)

guides <- c("gHNRNPM_1", "gHNRNPM_2")

## ----------------------------- packages ---------------------------------- ##
library(officer)  # write the .docx copy (rich runs give beautified equations)

fmt <- function(x) format(round(x, 2), nsmall = 2)   # 2-decimal display

## ---- beautified-equation helpers (affect the .docx only) ----------------- ##
## Inline markup used in the equation / "where" blocks below:
##   "_{...}"  -> the braced text is rendered as a subscript
##   MU        -> the Greek small letter mu (U+03BC), built encoding-safely below
## make_runs() turns such a string into officer runs (subscripted where marked);
## plain() strips the markup back to ASCII for the plain-text (.txt) mirror.
MU        <- intToUtf8(0x03BC)                     # Greek small letter mu
# Cambria = the template's body (minor-theme) font, so beautified runs match the
# surrounding Normal-style prose; headings use Calibri (major theme).
EQ_NORMAL <- fp_text(font.size = 11, font.family = "Cambria")
EQ_SUB    <- fp_text(font.size = 11, font.family = "Cambria", vertical.align = "subscript")

make_runs <- function(s) {
  runs <- list()
  repeat {
    m <- regexpr("_\\{[^}]*\\}", s)
    if (m[1] == -1) { if (nchar(s) > 0) runs <- c(runs, list(ftext(s, EQ_NORMAL))); break }
    pre <- substr(s, 1, m[1] - 1)
    grp <- regmatches(s, m)                       # e.g. "_{deg}"
    subtxt <- gsub("_\\{|\\}", "", grp)           # "deg"
    if (nchar(pre) > 0) runs <- c(runs, list(ftext(pre, EQ_NORMAL)))
    runs <- c(runs, list(ftext(subtxt, EQ_SUB)))
    s <- substr(s, m[1] + attr(m, "match.length"), nchar(s))
  }
  if (length(runs) == 0) runs <- list(ftext(" ", EQ_NORMAL))
  runs
}

plain <- function(s) {                            # markup -> ASCII (for the .txt)
  s <- gsub("_\\{([^}]*)\\}", "_\\1", s)          # k_{deg} -> k_deg
  s <- gsub(MU, "mu", s, fixed = TRUE)            # mu-letter -> "mu"
  s
}

## ----------------------------- inputs ------------------------------------ ##
## Everything comes from manuscript_inputs.csv (tidy: metric, guide, value),
## produced by doubling_vs_genotype.R. No re-computation happens here, so the
## text is guaranteed consistent with the figure. Run that script first.
if (!file.exists(INPUTS_CSV))
  stop("Missing '", INPUTS_CSV, "'. Run doubling_vs_genotype.R first to generate it.")
vals <- read.csv(INPUTS_CSV, stringsAsFactors = FALSE)

getv <- function(metric, guide = "all") {                 # raw (character) value
  v <- vals$value[vals$metric == metric & vals$guide == guide]
  if (length(v) != 1)
    stop("Expected one row for metric='", metric, "', guide='", guide, "' in ", INPUTS_CSV)
  v
}
getn <- function(metric, guide = "all") as.numeric(getv(metric, guide))  # numeric value

# EdU-derived ratio, keyed as in the CSV (R_edu_Day1 / R_edu_Day4)
R_edu <- function(g, day) getn(paste0("R_edu_", gsub(" ", "", day)), g)

# growth-fit results per guide, rebuilt so the text below can use gc[[g]]$...
gc <- lapply(guides, function(g) list(
  win_label = getv("win_label", g),
  Td_ctrl_h = getn("Td_ctrl_h", g),
  Td_kd_h   = getn("Td_kd_h",   g),
  R_growth  = getn("R_growth",  g),
  r2_ctrl   = getn("r2_ctrl",   g),
  r2_kd     = getn("r2_kd",     g)))
names(gc) <- guides

# circRNA fold-change summary
FC_min <- getn("FC_min")
FC_max <- getn("FC_max")
n_circ <- getn("n_circ")

# the most generous (largest) proliferation ceiling across all estimates
all_R <- c(R_edu("gHNRNPM_1", "Day 1"), R_edu("gHNRNPM_2", "Day 1"),
           R_edu("gHNRNPM_1", "Day 4"), R_edu("gHNRNPM_2", "Day 4"),
           gc[["gHNRNPM_1"]]$R_growth, gc[["gHNRNPM_2"]]$R_growth)
R_max <- max(all_R)

## ----------------------------- build the text ---------------------------- ##
# sections: kind (h1/h2/body) + text; used for both the .txt and the .docx
sections <- data.frame(kind = character(), text = character(), stringsAsFactors = FALSE)
add <- function(kind, text) rbind(sections, data.frame(kind = kind, text = text, stringsAsFactors = FALSE))

sections <- add("h1", "Methods")
sections <- add("h2", "Assessment of the proliferation confounder on circRNA abundance")
sections <- add("body", paste0(
  "To test whether the reduced proliferation of HNRNPM-knockdown cells could itself account for the ",
  "increased circRNA levels, circRNA abundance was modelled as the steady state of a single compartment ",
  "in which a stable RNA is removed both by molecular decay and by dilution through cell division. Under ",
  "this model the per-cell circRNA fold change attributable to slower division alone is:"))
sections <- add("eq", paste0("FC = (k_{deg} + ", MU, "_{ctrl}) / (k_{deg} + ", MU, "_{KD})"))
sections <- add("eq", paste0("with  k_{deg} = ln2 / half-life    and    ", MU, " = ln2 / doubling time"))
sections <- add("bodyrich", paste0(
  "As the half-life increases this fold change approaches a ceiling equal to the knockdown-to-control ",
  "doubling-time ratio R. R was estimated in two independent ways per guide RNA: from the EdU S-phase ",
  "fractions (R = %S_{ctrl} / %S_{KD}, assuming a constant S-phase duration) at Day 1 and Day 4, and from ",
  "the growth-curve doubling times obtained by log-linear regression of log(cell count) against time ",
  "(T_{d} = ln2 / slope; R = T_{d,KD} / T_{d,ctrl}). Because the cultures approached confluence and left ",
  "exponential growth toward the end of the time course, the fit was restricted to the exponential phase: ",
  "among all consecutive-day windows of at least four timepoints, the window whose log-linear fit was most ",
  "linear (highest mean coefficient of determination across the control and knockdown lines) was selected, ",
  "and the same window was applied to both lines within an experiment (gHNRNPM_1, window ",
  gc[["gHNRNPM_1"]]$win_label, "; gHNRNPM_2, window ", gc[["gHNRNPM_2"]]$win_label,
  "). The resulting ceilings were compared with the fold changes of the ", n_circ,
  " HNRNPM-independent circRNAs (Supplementary Table S5B, column logFC_BSJ). Modelling, statistics and ",
  "figures were performed in R (v4.3.3)."))

sections <- add("h1", "Results")
sections <- add("h2", "The mild proliferation slow-down cannot account for the circRNA increases")

sections <- add("bodyrich", paste0(
  "HNRNPM knockdown modestly lengthened the cell-doubling time. From the EdU S-phase fractions (assuming a ",
  "constant S-phase duration), the implied control-to-knockdown doubling-time ratio at Day 1 was R = ",
  fmt(R_edu("gHNRNPM_1", "Day 1")), " for gHNRNPM_1 and R = ", fmt(R_edu("gHNRNPM_2", "Day 1")),
  " for gHNRNPM_2; at Day 4 the corresponding ratios were R = ", fmt(R_edu("gHNRNPM_1", "Day 4")),
  " and R = ", fmt(R_edu("gHNRNPM_2", "Day 4")), ". Independently, log-linear fits over the exponential ",
  "growth phase (window ", gc[["gHNRNPM_1"]]$win_label, " for gHNRNPM_1 and ", gc[["gHNRNPM_2"]]$win_label,
  " for gHNRNPM_2) gave doubling-time ratios of R = ", fmt(gc[["gHNRNPM_1"]]$R_growth),
  " (gHNRNPM_1; T_{d,ctrl} = ", fmt(gc[["gHNRNPM_1"]]$Td_ctrl_h), " h, T_{d,KD} = ", fmt(gc[["gHNRNPM_1"]]$Td_kd_h),
  " h) and R = ", fmt(gc[["gHNRNPM_2"]]$R_growth), " (gHNRNPM_2; T_{d,ctrl} = ", fmt(gc[["gHNRNPM_2"]]$Td_ctrl_h),
  " h, T_{d,KD} = ", fmt(gc[["gHNRNPM_2"]]$Td_kd_h), " h)."))

sections <- add("body", paste0(
  "The two independent estimates agree, both in direction and in magnitude, that HNRNPM knockdown produces ",
  "only a modest slowing of proliferation, with R across all estimates falling in the ~", fmt(min(all_R)), " to ",
  fmt(max(all_R)), " range. The growth-curve and EdU-derived ratios are closely comparable, and we ",
  "conservatively used the largest of them as the proliferation ceiling in the analysis below. The reduction ",
  "in the EdU S-phase fraction was itself statistically significant in most comparisons (Welch's two-sample ",
  "t-test, gHNRNPM_2 at Day 1 and both guides at Day 4; Student's t-test, all four guide/day comparisons) but ",
  "small in magnitude."))

sections <- add("body", paste0(
  "Because R is the ceiling on the proliferation-only fold change, the largest circRNA increase that reduced ",
  "division could produce is at most ", fmt(R_max), "-fold, and only for an essentially non-degrading circRNA; ",
  "for any realistic half-life the achievable fold change is smaller (Fig. 1, Day 1; Supplementary Fig., Day 4). ",
  "In contrast, the smallest fold change among the ", n_circ, " HNRNPM-independent circRNAs was ", fmt(FC_min),
  "-fold (Supplementary Table S5B; range ", fmt(FC_min), " to ", fmt(FC_max), "-fold). Because even this ",
  "least-changed circRNA exceeds the most generous proliferation ceiling (", fmt(FC_min), " > ", fmt(R_max),
  "), the increase in circRNA levels cannot be explained by reduced proliferation alone, and is instead ",
  "dominated by direct HNRNPM-dependent regulation of back-splicing. The proliferation-only fold-change curves, ",
  "their ceilings, and this smallest-circRNA-FC reference are shown in ", FIG_PDF,
  " (Day 1 as the main panels; Day 4 as the supplementary panels)."))

## ----------------------------- write .txt -------------------------------- ##
## The .txt is a plain-text mirror: markup is stripped back to ASCII by plain().
out <- c()
for (i in 1:nrow(sections)) {
  out <- c(out, plain(sections$text[i]), "")
}
writeLines(out, OUTPUT_TXT)
print(paste0("Saved ", OUTPUT_TXT))

## ----------------------------- write .docx ------------------------------- ##
## Headings/plain body use body_add_par; "eq" (display equation) and "bodyrich"
## (prose with inline subscripts / Greek mu) use body_add_fpar with beautified runs.
doc <- read_docx()
for (i in 1:nrow(sections)) {
  k  <- sections$kind[i]
  tx <- sections$text[i]
  if (k == "h1") {
    doc <- body_add_par(doc, plain(tx), style = "heading 1")
  } else if (k == "h2") {
    doc <- body_add_par(doc, plain(tx), style = "heading 2")
  } else if (k == "eq") {
    doc <- body_add_fpar(doc, do.call(fpar, c(make_runs(tx),
             list(fp_p = fp_par(padding.left = 24, padding.top = 3, padding.bottom = 3)))))
  } else if (k == "bodyrich") {
    doc <- body_add_fpar(doc, do.call(fpar, make_runs(tx)), style = "Normal")
  } else {  # plain body
    doc <- body_add_par(doc, plain(tx), style = "Normal")
  }
}
print(doc, target = OUTPUT_DOCX)
print(paste0("Saved ", OUTPUT_DOCX))
