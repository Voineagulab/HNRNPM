# Fig S4A — circRNA x differential-AS overlap by AS class: BSJ-additive-DE (56) vs universe.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d1_circDE_vs_AS/d1_circDE_vs_AS.R
#   (the section that generates grouped_pct_overlap_56.pdf). Reads the two saved result
#   tables that d1_circDE_vs_AS.R wrote (Fisher-by-class and the per-stratum % summary)
#   rather than recomputing the overlap, and reproduces the grouped bar plot: for each
#   rMATS AS class, the % of circRNAs that overlap >=1 significant HNRNPM-differential AS
#   event, in the 56 additive BSJ-DE circRNAs (red) versus the whole detected universe (grey).
#
# Strategy (same as the Fig_2/Fig_3/Fig_4/Fig_5 series): keep ALL original plot parameters
#   (bar width, dodge, data-label sizes, OR/P annotation sizes, title/subtitle, bottom legend),
#   render at the original size, and let ggsave(scale=) squash the figure to 16 x 6.83 cm
#   (uniform shrink preserves proportions). Pre-inflate ONLY the axis title/tick fonts + legend
#   so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger than
#   nominal; resize it to 16 cm wide on the slide. See Fig_2E_exactsize.R.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
D1     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d1_circDE_vs_AS"
FISHER <- file.path(D1, "circRNA36_AS_overlap_fisher_byclass.tsv")
SUMM   <- file.path(D1, "circRNA36_AS_overlap_summary_byclass.tsv")
OUT_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/Figures/Supplementary"

SCOPE   <- "sig_full"                                  # significant AS events only
COMPARE <- "BSJ_DE vs universe"                        # Fisher comparison row (subset vs universe)
STRATUM <- "BSJ_DE"                                    # summary stratum for the 56 BSJ-DE circRNAs
ORDER   <- c("SE","MXE","A3SS","A5SS","RI","any")
LAB     <- function(c) ifelse(c == "any", "Any", c)
SIG <- "#d73027"; UNI <- "#4d4d4d"                     # red = BSJ-DE (56); grey = universe
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 16; H_CM <- 6.83          # nominal final PDF size on the slide
ORIG_W_IN <- 8                       # the source rendered grouped_pct_overlap_56.pdf at width = 8 in
SC     <- (ORIG_W_IN * 2.54) / W_CM     # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                         # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC       # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- read saved result tables ----
FI <- read.table(FISHER, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
SU <- read.table(SUMM,   sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

fi <- FI[FI$comparison == COMPARE & FI$scope == SCOPE, ]; rownames(fi) <- fi$AS_class
getp <- function(st, cl) { r <- SU[SU$scope == SCOPE & SU$stratum == st & SU$AS_class == cl, ]
  c(pct = r$pct_match[1], k = r$n_match[1], n = r$n_total[1]) }

D <- do.call(rbind, lapply(ORDER, function(cl) {
  ip <- getp(STRATUM, cl); up <- getp("universe", cl)
  data.frame(cl = cl, OR = fi[cl,"OR"], P = fi[cl,"p_one_sided"],
             sub_pct = ip["pct"], sub_k = ip["k"], sub_n = ip["n"],
             uni_pct = up["pct"], uni_k = up["k"], uni_n = up["n"], stringsAsFactors = FALSE) }))
D$cl <- factor(D$cl, levels = ORDER); D$xlab <- factor(LAB(as.character(D$cl)), levels = LAB(ORDER))
D$sig <- ifelse(D$P < 0.05, "sig", "ns")

red_label <- "BSJ-additive-DE (n=56)"
long <- rbind(
  data.frame(cl = D$cl, xlab = D$xlab, grp = red_label,  pct = D$sub_pct, k = D$sub_k, n = D$sub_n),
  data.frame(cl = D$cl, xlab = D$xlab, grp = "universe", pct = D$uni_pct, k = D$uni_k, n = D$uni_n))
long$grp <- factor(long$grp, levels = c(red_label, "universe"))
ymax <- max(long$pct)
ann  <- data.frame(xlab = D$xlab, y = ymax * 1.28, lab = sprintf("OR=%.2f\nP=%s", D$OR, fmt_p(D$P)), sig = D$sig)
fills <- c(SIG, UNI); names(fills) <- c(red_label, "universe")

# one continuous subtitle string, then wrapped at ~full-width (in characters) so each
# line fills the plot width before breaking to the next line (ggplot does not auto-wrap)
subtitle_txt <- paste0(
  "Per AS class: one-sided Fisher's exact test of overlap, BSJ-additive-DE vs universe; AS events significant at FDR<0.05 & |dPSI|>=0.10. ",
  "Bars = % of circRNAs from either the 56 BSJ-additive-DE group (red) or the whole detected set (grey) that overlap >=1 HNRNPM-regulated AS event. ",
  "k/n above bars: k = circRNAs overlapping, n = circRNAs in the group. OR = sample odds ratio.")
subtitle_txt <- paste(strwrap(subtitle_txt, width = 190), collapse = "\n")

# ---- build the grouped bar plot (original parameters preserved) ----
p <- ggplot(long, aes(xlab, pct, fill = grp)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d/%d", k, n)), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.4) +
  geom_text(data = ann, aes(x = xlab, y = y, label = lab), inherit.aes = FALSE, size = 2.6,
            colour = ifelse(ann$sig == "sig", "black", "grey55"), lineheight = 0.9) +
  scale_fill_manual(values = fills, name = NULL) +
  scale_y_continuous(limits = c(0, ymax * 1.45), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "AS event class", y = "% of circRNAs sharing splice sites with \n >=1 HNRNPM-regulated AS event",
       title = "circRNA x differential-AS overlap by class: BSJ-additive-DE (56) vs universe",
       subtitle = subtitle_txt) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title    = element_text(size = 5, margin = margin(b = 0)),   # font 5; no gap below title
        plot.subtitle = element_text(size = 5, margin = margin(t = 0)),   # font 5; no gap above subtitle
        axis.title    = element_text(size = axis_render_pt),  # -> 8 pt after shrink
        axis.text     = element_text(size = axis_render_pt),  # -> 8 pt after shrink
        legend.text   = element_text(size = axis_render_pt))  # -> 8 pt after shrink

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_S4A.pdf"), p,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_S4A.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_S4A.pdf)
# ============================================================
# Co-localisation of HNRNPM-BSJ-DE circRNAs with HNRNPM-regulated alternative splicing,
# broken down by alternative-splicing (AS) event class. For each rMATS event class
# (SE skipped exon, MXE mutually exclusive exons, A3SS/A5SS alternative 3'/5' splice site,
# RI retained intron, and "Any" = the union across all five classes), bars show the
# percentage of circRNAs whose back-splice splice site coincides with the splice site of
# at least one significant HNRNPM-differential AS event at the same host gene. Two groups
# are compared per class: the 56 additive (~ group + PSD) BSJ-differentially-expressed
# circRNAs (red) and the whole detected circRNA universe of 5,820 (grey). An AS event is
# called significant at rMATS FDR < 0.05 and |dPSI| >= 0.10. The count above each bar is
# k/n, where k is the number of circRNAs in that group overlapping such an event and n is
# the number of circRNAs in the group. Above each class the sample odds ratio (OR) and the
# one-sided Fisher's exact test P value (alternative, BSJ-DE more enriched than universe)
# are annotated, printed in grey where P >= 0.05. The horizontal axis is the AS event class
# and the vertical axis is the percentage of circRNAs overlapping a significant AS event.
# BSJ-DE circRNAs are enriched for co-localising with HNRNPM-differential AS events overall
# and specifically at skipped-exon (SE) and mutually-exclusive-exon (MXE) events, linking
# HNRNPM-dependent back-splicing to HNRNPM-dependent linear alternative splicing at shared
# splice sites. P values are per class and are not corrected for multiple testing.
