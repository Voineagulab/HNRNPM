#!/usr/bin/env Rscript
# =============================================================================
# d1_circDE_vs_AS_any.R
# -----------------------------------------------------------------------------
# ALL-EVENTS ("_ANY") counterpart of the two grouped-% panels in
#   d1_circDE_vs_AS.R  (same folder).
#
# WHAT DIFFERS FROM d1_circDE_vs_AS.R
#   d1_circDE_vs_AS.R plots SCOPE = "sig_full": circRNA splice sites matched only
#   against HNRNPM-regulated AS events (FDR < 0.05 & |dPSI| >= 0.10; 2,617 events).
#   This script plots SCOPE = "all": the same circRNA groups matched against
#   EVERY AS event detected by rMATS, with no FDR / dPSI filter (135,196 events
#   in .../FinalAnalysis/RESULTS/07_Splicing_rMATS/rmats_raw/).
#
# WHY NO NEW STATISTICS ARE COMPUTED HERE
#   The "all" scope was already computed upstream and is already stored:
#     - FinalAnalysis item 08 (circRNA_AS_junction_overlap_PSD58.py) read the raw
#       rmats_raw/*.MATS.JCEC.txt files and built BOTH a no-filter index (idx_all)
#       and a significant-only index (idx_sig_full), writing matches_*_all and
#       matches_*_sig_full flags into circRNA_AS_overlap_PSD58.tsv.
#     - d1_circDE_vs_AS.R then looped over SCOPES <- c("all", "sig_full") for the
#       Fisher tests and the per-stratum % summary, so both of its TSVs already
#       carry 'all' rows (18 of 36 in the Fisher table, 24 of 48 in the summary).
#   This script therefore READS those saved rows rather than recomputing them, and
#   writes NO new TSV. Nothing outside this folder is read or written.
#
# READ THE PANELS AS A CONTROL, NOT AS A RESULT
#   The "all" baseline is saturated: 95.6% of the 5,820-circRNA universe already
#   sits at the splice site of >=1 rMATS-tested event, as do 36/36 independent-36
#   and 55/56 BSJ-DE circRNAs. There is thus no headroom for an enrichment test.
#   Every one-sided P in the 'all' scope is >= 0.095, and the independent-36 odds
#   ratio is Inf for SE and for 'any' because its non-matching cell is zero (that
#   label renders literally as "OR=Inf", the same sprintf("%.2f") behaviour as the
#   source script -- left unchanged deliberately). These panels are useful only to
#   show that the enrichment seen at SCOPE = "sig_full" is specific to
#   HNRNPM-regulated events.
#
# OUTPUTS (this folder only)
#   grouped_pct_overlap_ANY.pdf     independent-36 vs universe, all detected events
#   grouped_pct_overlap_56_ANY.pdf  BSJ-additive-DE (56) vs universe, all events
#
# BAR FILL: the _ANY bars use a LIGHT TINT of the same two colours rather than the
#   full-strength red / dark grey (see tint() below). Pattern fills were tried first
#   and silently did not render - see the note above tint().
# PLOT HEIGHT: every ggsave() height in this script is 70% of its earlier value.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

D1   <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d1_circDE_vs_AS"
FI_F <- file.path(D1, "circRNA36_AS_overlap_fisher_byclass.tsv")    # written by d1_circDE_vs_AS.R
SU_F <- file.path(D1, "circRNA36_AS_overlap_summary_byclass.tsv")   # written by d1_circDE_vs_AS.R

for (f in c(FI_F, SU_F)) if (!file.exists(f))
  stop("Required input not found (run d1_circDE_vs_AS.R first): ", f)

FI <- read.table(FI_F, sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                 stringsAsFactors = FALSE)
SU <- read.table(SU_F, sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                 stringsAsFactors = FALSE)

# ---- the one substantive change vs d1_circDE_vs_AS.R ----
SCOPE <- "all"                       # d1_circDE_vs_AS.R uses "sig_full"

if (!any(FI$scope == SCOPE)) stop("No rows with scope == '", SCOPE, "' in ", FI_F)
if (!any(SU$scope == SCOPE)) stop("No rows with scope == '", SCOPE, "' in ", SU_F)
cat(sprintf("scope = '%s': %d Fisher rows, %d summary rows\n",
            SCOPE, sum(FI$scope == SCOPE), sum(SU$scope == SCOPE)))

# ---- everything below is the plotting code of d1_circDE_vs_AS.R, unchanged ----
ORDER <- c("SE", "MXE", "A3SS", "A5SS", "RI", "any")
LAB   <- function(c) ifelse(c == "any", "Any", c)
SIG <- "#d73027"; NS <- "#bdbdbd"; UNI <- "#4d4d4d"
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1),
                            formatC(p, format = "f", digits = 3))

# ---- distinguishing the two scopes: a light tint of the same two colours ----
# Earlier versions of this script tried grid pattern fills (diagonal hatch, then
# dots) supplied through scale_fill_manual(). That FAILED SILENTLY: ggplot2 3.5.1
# took the colour and discarded the pattern, so the bars rendered solid and the two
# scopes were indistinguishable while the legend still claimed a difference. Verified
# by sampling the 600 dpi raster - zero white pixels inside any bar. Patterns would
# need the ggpattern package, which is not installed here. Plain tinted colours
# cannot fail this way, so that is what is used.
#
# tint() moves a colour TINT of the way toward white: TINT = 0.55 leaves 45% of the
# original, giving #d73027 -> a pale red and #4d4d4d -> a mid grey. Raise TINT for
# paler ANY bars, lower it for a subtler distinction.
TINT <- 0.55

tint <- function(col, amount = TINT) {
  v <- grDevices::col2rgb(col) / 255
  v <- v + (1 - v) * amount
  grDevices::rgb(v[1], v[2], v[3])
}

make_grouped <- function(compare, stratum, red_label, title, subtitle, out_pdf) {
  fi <- FI[FI$comparison == compare & FI$scope == SCOPE, ]; rownames(fi) <- fi$AS_class
  if (!nrow(fi)) stop("No Fisher rows for comparison '", compare, "' at scope '", SCOPE, "'")
  getp <- function(st, cl) { r <- SU[SU$scope == SCOPE & SU$stratum == st & SU$AS_class == cl, ]
    c(pct = r$pct_match[1], k = r$n_match[1], n = r$n_total[1]) }
  D <- do.call(rbind, lapply(ORDER, function(cl) {
    ip <- getp(stratum, cl); up <- getp("universe", cl)
    data.frame(cl = cl, OR = fi[cl, "OR"], P = fi[cl, "p_one_sided"],
               sub_pct = ip["pct"], sub_k = ip["k"], sub_n = ip["n"],
               uni_pct = up["pct"], uni_k = up["k"], uni_n = up["n"], stringsAsFactors = FALSE) }))
  D$cl <- factor(D$cl, levels = ORDER); D$xlab <- factor(LAB(as.character(D$cl)), levels = LAB(ORDER))
  D$sig <- ifelse(D$P < 0.05, "sig", "ns")
  long <- rbind(
    data.frame(cl = D$cl, xlab = D$xlab, grp = red_label, pct = D$sub_pct, k = D$sub_k, n = D$sub_n),
    data.frame(cl = D$cl, xlab = D$xlab, grp = "universe", pct = D$uni_pct, k = D$uni_k, n = D$uni_n))
  long$grp <- factor(long$grp, levels = c(red_label, "universe"))
  ymax <- max(long$pct)
  ann <- data.frame(xlab = D$xlab, y = ymax * 1.28,
                    lab = sprintf("OR=%.2f\nP=%s", D$OR, fmt_p(D$P)), sig = D$sig)
  # Light-tinted rather than full-strength, keeping red for the circRNA group and
  # grey for the universe (these are the _ANY panels, i.e. all detected events).
  fills <- c(tint(SIG), tint(UNI)); names(fills) <- c(red_label, "universe")
  p <- ggplot(long, aes(xlab, pct, fill = grp)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d/%d", k, n)), position = position_dodge(width = 0.8),
              vjust = -0.3, size = 2.4) +
    geom_text(data = ann, aes(x = xlab, y = y, label = lab), inherit.aes = FALSE, size = 2.6,
              colour = ifelse(ann$sig == "sig", "black", "grey55"), lineheight = 0.9,
              fontface = "bold") +
    scale_fill_manual(values = fills, name = NULL) +
    scale_y_continuous(limits = c(0, ymax * 1.45), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "AS event class", y = "% of circRNAs overlapping a detected AS event",
         title = title, subtitle = subtitle) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 6))
  ggsave(out_pdf, p, width = 8, height = 3.64, device = cairo_pdf)   # height = 70% of the previous 5.2
  cat(sprintf("Wrote %s\n", out_pdf))
}

# independent-36 vs universe, ALL detected AS events
make_grouped(
  compare   = "independent36_DE vs universe", stratum = "independent36_DE",
  red_label = "independent-36",
  title     = "circRNA x ALL detected-AS overlap by class: independent-36 vs universe",
  subtitle  = paste0(
    "Per AS class: one-sided Fisher's exact test of overlap, independent-36 vs universe; ALL AS events detected by rMATS (no FDR / |dPSI| filter).\n",
    "Bars = % of circRNAs from either the 36-independent group (red) or the whole detected set (grey) that overlap >=1 detected AS event.\n",
    "k/n above bars: k = circRNAs overlapping, n = circRNAs in the group. OR = sample odds ratio."),
  out_pdf   = file.path(D1, "grouped_pct_overlap_ANY.pdf"))

# BSJ-additive-DE (56) vs universe, ALL detected AS events
make_grouped(
  compare   = "BSJ_DE vs universe", stratum = "BSJ_DE",   # "_additive" was stripped by d1_circDE_vs_AS.R
  red_label = "BSJ-additive-DE (n=56)",
  title     = "circRNA x ALL detected-AS overlap by class: BSJ-additive-DE (56) vs universe",
  subtitle  = paste0(
    "Per AS class: one-sided Fisher's exact test of overlap, BSJ-additive-DE vs universe; ALL AS events detected by rMATS (no FDR / |dPSI| filter).\n",
    "Bars = % of circRNAs from either the 56 BSJ-additive-DE group (red) or the whole detected set (grey) that overlap >=1 detected AS event.\n",
    "k/n above bars: k = circRNAs overlapping, n = circRNAs in the group. OR = sample odds ratio."),
  out_pdf   = file.path(D1, "grouped_pct_overlap_56_ANY.pdf"))

cat(sprintf("\nAll outputs written to %s\n", D1))

# =============================================================================
# MERGED PANELS - the two scopes side by side in one plot
# -----------------------------------------------------------------------------
# For each AS class, four bars: the pair from the sig_full panel (solid) next to
# the pair from the 'all' panel (light-tinted). Same two hues throughout
# (red = circRNA group, dark grey = universe); only the fill style distinguishes
# the scope.
#
# The scope distinction is carried by fill lightness, not by a pattern: see the
# note above tint() for why pattern fills were abandoned.
#
# Each class carries two OR / P annotations, "sig:" for the FDR<0.05 & |dPSI|>=0.10
# scope and "all:" for the unfiltered scope; the block is printed in black if
# either scope is significant, grey otherwise. OR renders as "Inf" where the
# group's non-matching cell is zero (independent-36 at SE and any, all scope) -
# the same sprintf("%.2f") behaviour as the source script.
#
# Wider than the single-scope panels (10 in vs 8 in) because there are now four
# bars per class instead of two.
#
# OUTPUTS
#   grouped_pct_overlap_ANY_merge.pdf     independent-36:    sig_full + all
#   grouped_pct_overlap_56_ANY_merge.pdf  BSJ-additive-DE 56: sig_full + all
# =============================================================================

SCOPE_SIG <- "sig_full"      # HNRNPM-regulated AS events (FDR<0.05 & |dPSI|>=0.10)
SCOPE_ANY <- "all"           # every AS event detected by rMATS

make_merged <- function(compare, stratum, red_label, title, subtitle, out_pdf) {
  getf <- function(sc, cl) {
    r <- FI[FI$scope == sc & FI$comparison == compare & FI$AS_class == cl, ]
    if (!nrow(r)) stop("No Fisher row: scope=", sc, " comparison=", compare, " class=", cl)
    c(OR = r$OR[1], P = r$p_one_sided[1])
  }
  getp <- function(sc, st, cl) {
    r <- SU[SU$scope == sc & SU$stratum == st & SU$AS_class == cl, ]
    if (!nrow(r)) stop("No summary row: scope=", sc, " stratum=", st, " class=", cl)
    c(pct = r$pct_match[1], k = r$n_match[1], n = r$n_total[1])
  }

  # Universe size taken from the summary table itself so it cannot go stale
  # (5,820 circRNAs detected in the SH-SY5Y samples).
  uni_n   <- unname(getp(SCOPE_SIG, "universe", "any")["n"])
  UNI_LAB <- sprintf("universe (n=%s)", format(uni_n, big.mark = ","))
  G <- c(paste0(red_label, " - HNRNPM-regulated AS"),
         paste0(UNI_LAB,   " - HNRNPM-regulated AS"),
         paste0(red_label, " - all detected AS"),
         paste0(UNI_LAB,   " - all detected AS"))

  long <- do.call(rbind, lapply(ORDER, function(cl) {
    v <- list(getp(SCOPE_SIG, stratum, cl), getp(SCOPE_SIG, "universe", cl),
              getp(SCOPE_ANY, stratum, cl), getp(SCOPE_ANY, "universe", cl))
    data.frame(cl  = cl, grp = G,
               pct = unname(sapply(v, function(x) x["pct"])),
               k   = unname(sapply(v, function(x) x["k"])),
               n   = unname(sapply(v, function(x) x["n"])),
               stringsAsFactors = FALSE) }))
  long$xlab <- factor(LAB(long$cl), levels = LAB(ORDER))
  long$grp  <- factor(long$grp, levels = G)
  # Only k above each bar; the group sizes (n) are given in the legend.
  long$kn   <- as.character(long$k)

  ymax <- max(long$pct, na.rm = TRUE)
  # One OR/P label per PAIR of bars, centred over that pair and lifted the SAME
  # fixed clearance (GAP) above the taller bar of the pair. So each label sits
  # directly above its own two bars, and the two labels in a class are no longer
  # forced onto one shared horizontal line.
  DODGE   <- 0.85                                   # must match position_dodge() below
  PAIR_DX <- c(sig = -DODGE / 4, any = DODGE / 4)    # pair centres, offset from class centre
  GAP     <- 0.22 * ymax                             # equal clearance for every pair (2x the previous 0.11)
  ann <- do.call(rbind, lapply(seq_along(ORDER), function(i) {
    cl <- ORDER[i]
    fs <- getf(SCOPE_SIG, cl); fa <- getf(SCOPE_ANY, cl)
    top_sig <- max(getp(SCOPE_SIG, stratum, cl)["pct"], getp(SCOPE_SIG, "universe", cl)["pct"])
    top_any <- max(getp(SCOPE_ANY, stratum, cl)["pct"], getp(SCOPE_ANY, "universe", cl)["pct"])
    data.frame(
      x   = i + unname(PAIR_DX[c("sig", "any")]),
      y   = unname(c(top_sig, top_any)) + GAP,
      lab = c(sprintf("OR=%.2f\nP=%s", fs["OR"], fmt_p(fs["P"])),
              sprintf("OR=%.2f\nP=%s", fa["OR"], fmt_p(fa["P"]))),
      sig = c(ifelse(fs["P"] < 0.05, "sig", "ns"), ifelse(fa["P"] < 0.05, "sig", "ns")),
      stringsAsFactors = FALSE) }))
  Y_TOP <- max(ann$y) + 0.16 * ymax                  # headroom for the two-line labels

  fills <- c(SIG, UNI, tint(SIG), tint(UNI)); names(fills) <- G

  p <- ggplot(long, aes(xlab, pct, fill = grp)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.78,
             colour = "black", linewidth = 0.3) +
    geom_text(aes(label = kn),
              position = position_dodge(width = 0.85), vjust = -0.3, size = 2.0) +
    geom_text(data = ann, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
              size = 2.4, colour = ifelse(ann$sig == "sig", "black", "grey55"),
              lineheight = 0.9, fontface = "bold") +
    scale_fill_manual(values = fills, name = NULL) +
    scale_y_continuous(limits = c(0, Y_TOP), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "AS event class", y = "% of circRNAs overlapping an AS event",
         title = title, subtitle = subtitle) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 6)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  ggsave(out_pdf, p, width = 10, height = 3.92, device = cairo_pdf)  # height = 70% of the previous 5.6
  cat(sprintf("Wrote %s\n", out_pdf))
}

MERGE_SUB <- paste0(
  "Full-strength bars = AS events significant at FDR<0.05 & |dPSI|>=0.10 (as in grouped_pct_overlap%s.pdf); pale bars = ALL AS events detected by rMATS (as in grouped_pct_overlap%s_ANY.pdf).\n",
  "Red = the circRNA group, grey = the whole detected set (universe). The number above each bar counts circRNAs overlapping >=1 event; group sizes (n) are in the legend.\n",
  "Per AS class, a one-sided Fisher's exact test vs universe is annotated for each scope; OR = sample odds ratio.\n",
  "NOTE: the unfiltered scope is saturated, so its ORs sit near 1 (or Inf where no non-matching circRNA remains) and are not evidence of enrichment.")

# independent-36: grouped_pct_overlap.pdf + grouped_pct_overlap_ANY.pdf
make_merged(
  compare   = "independent36_DE vs universe", stratum = "independent36_DE",
  red_label = "independent-36",
  title     = "circRNA x AS overlap by class, both scopes: independent-36 vs universe",
  subtitle  = sprintf(MERGE_SUB, "", ""),
  out_pdf   = file.path(D1, "grouped_pct_overlap_ANY_merge.pdf"))

# BSJ-additive-DE (56): grouped_pct_overlap_56.pdf + grouped_pct_overlap_56_ANY.pdf
make_merged(
  compare   = "BSJ_DE vs universe", stratum = "BSJ_DE",
  red_label = "BSJ-additive-DE (n=56)",
  title     = "circRNA x AS overlap by class, both scopes: BSJ-additive-DE (56) vs universe",
  subtitle  = sprintf(MERGE_SUB, "_56", "_56"),
  out_pdf   = file.path(D1, "grouped_pct_overlap_56_ANY_merge.pdf"))

cat(sprintf("\nMerged panels written to %s\n", D1))
