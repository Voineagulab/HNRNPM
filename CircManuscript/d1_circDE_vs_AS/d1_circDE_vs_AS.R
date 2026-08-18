#!/usr/bin/env Rscript
# =============================================================================
# d1_circDE_vs_AS.R
# -----------------------------------------------------------------------------
# Consolidated, self-contained copy of two scripts from
#   .../JuliWang_circRBP_pilot/z260721_36circDE_vs_AS_subclass_R/ :
#
#   PART 1  = full copy of circRNA36_AS_junction_overlap_byclass.R
#             -> writes circRNA36_AS_overlap_fisher_byclass.tsv
#                       circRNA36_AS_overlap_summary_byclass.tsv
#
#   PART 2  = partial copy of barplot/barplot_alt_representations.R
#             -> writes grouped_pct_overlap.pdf     (independent-36 vs universe)
#                       grouped_pct_overlap_56.pdf  (BSJ-additive-DE (56) vs universe)
#             (only the two grouped-% plots; the forest / stacked / combo panels
#              are NOT reproduced here)
#
# All outputs are written to this script's own folder
#   .../circRBP_pilot/IV/JW/CircManuscript/d1_circDE_vs_AS/
# Inputs are read (read-only) from the FinalAnalysis item-08 overlap table and the
# additive BSJ x FSJ classification. Nothing outside d1_circDE_vs_AS/ is modified.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

FA  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS"
OV  <- file.path(FA, "08_circRNA_AS_overlap", "circRNA_AS_overlap_PSD58.tsv")
CLS <- file.path(FA, "01_circRNA_DE", "FSJ", "BSJ_FSJ_classification_additive.tsv")
OUT <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d1_circDE_vs_AS"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EVENTS <- c("SE", "MXE", "A3SS", "A5SS", "RI")
SCOPES <- c("all", "sig_full")

ov  <- read.table(OV,  sep = "\t", header = TRUE, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
cls <- read.table(CLS, sep = "\t", header = TRUE, quote = "", check.names = FALSE, stringsAsFactors = FALSE)

indep <- cls$circRNA[cls$class == "independent_backsplicing"]
ov$INDEP_36        <- ov$circRNA %in% indep
ov$BSJ_additive_DE <- as.logical(ov$BSJ_additive_DE)
ov$FSJ_additive_DE <- as.logical(ov$FSJ_additive_DE)
cat(sprintf("independent_backsplicing circRNAs: %d (in universe: %d)\n", length(indep), sum(ov$INDEP_36)))

# =============================================================================
# PART 1 — Fisher by AS class (full copy of circRNA36_AS_junction_overlap_byclass.R)
# =============================================================================
strata <- list(independent36_DE = "INDEP_36",
               BSJ_additive_DE  = "BSJ_additive_DE",
               FSJ_additive_DE  = "FSJ_additive_DE")
matchcol <- function(cl, sc) if (cl == "any") paste0("matches_any_", sc) else paste0("matches_", cl, "_", sc)

fisher_one <- function(col, mask) {
  m  <- as.logical(ov[[col]])
  a  <- sum(m[mask]); sub_n <- sum(mask)
  M  <- sum(m); N <- length(m)
  b  <- sub_n - a
  cc <- M - a
  dd <- (N - sub_n) - cc
  tab <- matrix(c(a, b, cc, dd), nrow = 2, byrow = TRUE)
  p1 <- fisher.test(tab, alternative = "greater")$p.value
  p2 <- fisher.test(tab, alternative = "two.sided")$p.value
  OR <- (a * dd) / (b * cc)                                 # sample OR (matches scipy fisher_exact)
  list(sub_match = a, sub_n = sub_n, bg_match = cc, bg_n = N - sub_n, OR = OR, p1 = p1, p2 = p2)
}

fr <- list()
for (sc in SCOPES) for (cl in c(EVENTS, "any")) for (nm in names(strata)) {
  r <- fisher_one(matchcol(cl, sc), ov[[strata[[nm]]]])
  fr[[length(fr) + 1]] <- data.frame(scope = sc, AS_class = cl, comparison = paste(nm, "vs universe"),
    sub_match = r$sub_match, sub_n = r$sub_n, bg_match = r$bg_match, bg_n = r$bg_n,
    OR = round(r$OR, 3), p_one_sided = r$p1, p_two_sided = r$p2, stringsAsFactors = FALSE)
}
fisher <- do.call(rbind, fr)
# remove "_additive" from every string cell (e.g. BSJ_additive_DE -> BSJ_DE, FSJ_additive_DE -> FSJ_DE)
fisher[] <- lapply(fisher, function(x) if (is.character(x)) gsub("_additive", "", x, fixed = TRUE) else x)
write.table(fisher, file.path(OUT, "circRNA36_AS_overlap_fisher_byclass.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sr <- list(); uni <- rep(TRUE, nrow(ov))
smstrata <- list(independent36_DE = ov$INDEP_36,
                 BSJ_additive_DE  = ov$BSJ_additive_DE,
                 FSJ_additive_DE  = ov$FSJ_additive_DE,
                 universe         = uni)
for (sc in SCOPES) for (cl in c("any", EVENTS)) for (nm in names(smstrata)) {
  m <- as.logical(ov[[matchcol(cl, sc)]]); mask <- smstrata[[nm]]
  n <- sum(mask); k <- sum(m[mask])
  sr[[length(sr) + 1]] <- data.frame(stratum = nm, scope = sc, AS_class = cl,
    n_total = n, n_match = k, pct_match = round(100 * k / max(n, 1), 2), stringsAsFactors = FALSE)
}
summary <- do.call(rbind, sr)
# remove "_additive" from every string cell (e.g. BSJ_additive_DE -> BSJ_DE, FSJ_additive_DE -> FSJ_DE)
summary[] <- lapply(summary, function(x) if (is.character(x)) gsub("_additive", "", x, fixed = TRUE) else x)
write.table(summary, file.path(OUT, "circRNA36_AS_overlap_summary_byclass.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== independent36_DE vs universe (Fisher by AS class) ===\n")
print(fisher[fisher$comparison == "independent36_DE vs universe",
             c("scope","AS_class","sub_match","sub_n","OR","p_one_sided","p_two_sided")], row.names = FALSE)

# =============================================================================
# PART 2 — grouped % overlap plots (partial copy of barplot_alt_representations.R)
# =============================================================================
FI <- fisher; SU <- summary
SCOPE <- "sig_full"; ORDER <- c("SE","MXE","A3SS","A5SS","RI","any")
LAB <- function(c) ifelse(c == "any", "Any", c)
SIG <- "#d73027"; NS <- "#bdbdbd"; UNI <- "#4d4d4d"
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

make_grouped <- function(compare, stratum, red_label, title, subtitle, out_pdf) {
  fi <- FI[FI$comparison == compare & FI$scope == SCOPE, ]; rownames(fi) <- fi$AS_class
  getp <- function(st, cl) { r <- SU[SU$scope==SCOPE & SU$stratum==st & SU$AS_class==cl, ]
    c(pct = r$pct_match[1], k = r$n_match[1], n = r$n_total[1]) }
  D <- do.call(rbind, lapply(ORDER, function(cl) {
    ip <- getp(stratum, cl); up <- getp("universe", cl)
    data.frame(cl = cl, OR = fi[cl,"OR"], P = fi[cl,"p_one_sided"],
               sub_pct = ip["pct"], sub_k = ip["k"], sub_n = ip["n"],
               uni_pct = up["pct"], uni_k = up["k"], uni_n = up["n"], stringsAsFactors = FALSE) }))
  D$cl <- factor(D$cl, levels = ORDER); D$xlab <- factor(LAB(as.character(D$cl)), levels = LAB(ORDER))
  D$sig <- ifelse(D$P < 0.05, "sig", "ns")
  long <- rbind(
    data.frame(cl = D$cl, xlab = D$xlab, grp = red_label, pct = D$sub_pct, k = D$sub_k, n = D$sub_n),
    data.frame(cl = D$cl, xlab = D$xlab, grp = "universe", pct = D$uni_pct, k = D$uni_k, n = D$uni_n))
  long$grp <- factor(long$grp, levels = c(red_label, "universe"))
  ymax <- max(long$pct)
  ann <- data.frame(xlab = D$xlab, y = ymax * 1.28, lab = sprintf("OR=%.2f\nP=%s", D$OR, fmt_p(D$P)), sig = D$sig)
  fills <- c(SIG, UNI); names(fills) <- c(red_label, "universe")
  p <- ggplot(long, aes(xlab, pct, fill = grp)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d/%d", k, n)), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.4) +
    geom_text(data = ann, aes(x = xlab, y = y, label = lab), inherit.aes = FALSE, size = 2.6,
              colour = ifelse(ann$sig == "sig", "black", "grey55"), lineheight = 0.9) +
    scale_fill_manual(values = fills, name = NULL) +
    scale_y_continuous(limits = c(0, ymax * 1.45), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "AS event class", y = "% of circRNAs overlapping a significant AS event",
         title = title, subtitle = subtitle) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 6))
  ggsave(out_pdf, p, width = 8, height = 3.64, device = cairo_pdf)   # height = 70% of the original 5.2
  cat(sprintf("Wrote %s\n", out_pdf))
}

# independent-36 vs universe
make_grouped(
  compare   = "independent36_DE vs universe", stratum = "independent36_DE",
  red_label = "independent-36",
  title     = "circRNA x differential-AS overlap by class: independent-36 vs universe",
  subtitle  = paste0(
    "Per AS class: one-sided Fisher's exact test of overlap, independent-36 vs universe; AS events significant at FDR<0.05 & |dPSI|>=0.10.\n",
    "Bars = % of circRNAs from either the 36-independent group (red) or the whole detected set (grey) that overlap >=1 HNRNPM-regulated AS event.\n",
    "k/n above bars: k = circRNAs overlapping, n = circRNAs in the group. OR = sample odds ratio."),
  out_pdf   = file.path(OUT, "grouped_pct_overlap.pdf"))

# BSJ-additive-DE (56) vs universe
make_grouped(
  compare   = "BSJ_DE vs universe", stratum = "BSJ_DE",   # both fisher comparison and summary stratum had "_additive" stripped
  red_label = "BSJ-additive-DE (n=56)",
  title     = "circRNA x differential-AS overlap by class: BSJ-additive-DE (56) vs universe",
  subtitle  = paste0(
    "Per AS class: one-sided Fisher's exact test of overlap, BSJ-additive-DE vs universe; AS events significant at FDR<0.05 & |dPSI|>=0.10.\n",
    "Bars = % of circRNAs from either the 56 BSJ-additive-DE group (red) or the whole detected set (grey) that overlap >=1 HNRNPM-regulated AS event.\n",
    "k/n above bars: k = circRNAs overlapping, n = circRNAs in the group. OR = sample odds ratio."),
  out_pdf   = file.path(OUT, "grouped_pct_overlap_56.pdf"))

cat(sprintf("\nAll outputs written to %s\n", OUT))
