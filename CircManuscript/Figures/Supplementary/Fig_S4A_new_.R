# Fig S4A (new) — circRNA x AS overlap by class, BOTH scopes: BSJ-additive-DE (56) vs universe.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/Figures/Supplementary/Fig_S4A.R
#   (structure and house conventions) and from
#   /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d1_circDE_vs_AS/d1_circDE_vs_AS_any.R
#   (the make_merged() section that generates grouped_pct_overlap_56_ANY_merge.pdf, which this
#   script reproduces exactly). Reads the two saved result tables that d1_circDE_vs_AS.R wrote
#   (Fisher-by-class and the per-stratum % summary) rather than recomputing the overlap, and
#   draws FOUR bars per rMATS AS class: the significant-AS pair at full strength beside the
#   all-detected-AS pair in a light tint, for the 56 additive BSJ-DE circRNAs (red) and the
#   whole detected universe (grey).
#
# Strategy (same as the Fig_2/Fig_3/Fig_4/Fig_5 series): keep ALL original plot parameters
#   (bar width, dodge, data-label sizes, OR/P annotation sizes, title/subtitle, bottom legend),
#   render at the original size, and let ggsave(scale=) squash the figure to 16 x 6.272 cm
#   (uniform shrink preserves proportions). H_CM is DERIVED from the source panel's 10 x 3.92 in
#   aspect ratio, so the raw PDF has exactly the geometry of grouped_pct_overlap_56_ANY_merge.pdf.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger than
#   nominal; resize it to 16 cm wide on the slide. See Fig_2E_exactsize.R.
#
# DEPARTURES from the source panel grouped_pct_overlap_56_ANY_merge.pdf, which this
# script was originally written to reproduce exactly. It deliberately no longer does:
#   1. FONT SIZES follow the house convention rather than the source panel. legend.text,
#      axis.title.x and axis.text.x render at AXIS_PT * SC = 12.7 pt so they land at 8 pt once
#      the PDF is placed at 16 x 6.272 cm on the slide. NOTE the y-axis was left at the
#      theme_classic default by request, so its title and tick labels land at about 6.9 and
#      5.5 pt; change axis.title.x / axis.text.x to axis.title / axis.text below to bring
#      both axes to 8 pt.
#   2. plot.title is set to the subtitle's 6 pt rather than the theme_classic default.
#   3. The y-axis label reads "% of circRNAs sharing / splice sites with AS event".
#   4. The subtitle keeps the source panel's EXPLICIT line breaks rather than the
#      strwrap(width = 190) idiom, which would re-wrap the text and change the render.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
D1     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/d1_circDE_vs_AS"
FISHER <- file.path(D1, "circRNA36_AS_overlap_fisher_byclass.tsv")
SUMM   <- file.path(D1, "circRNA36_AS_overlap_summary_byclass.tsv")
OUT_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/Figures/Supplementary"

SCOPE_SIG <- "sig_full"                                # HNRNPM-regulated AS events only
SCOPE_ANY <- "all"                                     # every AS event detected by rMATS
COMPARE <- "BSJ_DE vs universe"                        # Fisher comparison row (subset vs universe)
STRATUM <- "BSJ_DE"                                    # summary stratum for the 56 BSJ-DE circRNAs
ORDER   <- c("SE","MXE","A3SS","A5SS","RI","any")
LAB     <- function(c) ifelse(c == "any", "Any", c)
SIG <- "#d73027"; UNI <- "#4d4d4d"                     # red = BSJ-DE (56); grey = universe
TINT <- 0.55                                           # fraction toward white for the all-scope bars
tint <- function(col, amount = TINT) {                 # #d73027 -> pale red; #4d4d4d -> mid grey
  v <- grDevices::col2rgb(col) / 255
  v <- v + (1 - v) * amount
  grDevices::rgb(v[1], v[2], v[3])
}
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM      <- 16                      # nominal final PDF width on the slide
ORIG_W_IN <- 10                      # the source rendered grouped_pct_overlap_56_ANY_merge.pdf at 10 in
ORIG_H_IN <- 3.92                    # ... and at 3.92 in high
H_CM      <- W_CM * ORIG_H_IN / ORIG_W_IN   # 6.272 cm; preserves the source aspect ratio exactly
SC        <- (ORIG_W_IN * 2.54) / W_CM      # ggsave scale: draw at original size, shrink to target
AXIS_PT   <- 8                       # desired FINAL legend + x-axis font size at W_CM x H_CM
axis_render_pt <- AXIS_PT * SC       # = 12.7 pt; pre-inflate so it lands at 8 pt after the shrink

# ---- read saved result tables ----
FI <- read.table(FISHER, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
SU <- read.table(SUMM,   sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

getf <- function(sc, cl) { r <- FI[FI$scope == sc & FI$comparison == COMPARE & FI$AS_class == cl, ]
  if (!nrow(r)) stop("No Fisher row: scope=", sc, " class=", cl)
  c(OR = r$OR[1], P = r$p_one_sided[1]) }
getp <- function(sc, st, cl) { r <- SU[SU$scope == sc & SU$stratum == st & SU$AS_class == cl, ]
  if (!nrow(r)) stop("No summary row: scope=", sc, " stratum=", st, " class=", cl)
  c(pct = r$pct_match[1], k = r$n_match[1], n = r$n_total[1]) }

# Universe size read from the summary table itself so it cannot go stale (5,820 circRNAs).
uni_n     <- unname(getp(SCOPE_SIG, "universe", "any")["n"])
sub_n     <- unname(getp(SCOPE_SIG, STRATUM,    "any")["n"])   # 56 BSJ-additive-DE circRNAs
UNI_LAB   <- sprintf("universe (n=%s)", format(uni_n, big.mark = ","))
red_label <- "BSJ-additive-DE (n=56)"
G <- c("DE circRNAs overlapping HNRNPM-regulated AS",
       "All circRNAs overlapping HNRNPM-regulated AS",
       "DE circRNAs overlapping all detected AS",
       "All circRNAs overlapping all detected AS")

long <- do.call(rbind, lapply(ORDER, function(cl) {
  v <- list(getp(SCOPE_SIG, STRATUM, cl), getp(SCOPE_SIG, "universe", cl),
            getp(SCOPE_ANY, STRATUM, cl), getp(SCOPE_ANY, "universe", cl))
  data.frame(cl  = cl, grp = G,
             pct = unname(sapply(v, function(x) x["pct"])),
             k   = unname(sapply(v, function(x) x["k"])),
             n   = unname(sapply(v, function(x) x["n"])),
             stringsAsFactors = FALSE) }))
long$xlab <- factor(LAB(long$cl), levels = LAB(ORDER))
long$grp  <- factor(long$grp, levels = G)
long$kn   <- as.character(long$k)      # only k above each bar; group sizes (n) are in the legend

ymax <- max(long$pct, na.rm = TRUE)
# One OR/P label per PAIR of bars, centred over that pair and lifted the SAME fixed clearance
# (GAP) above the taller bar of the pair, so the two labels in a class are not forced onto one
# shared horizontal line.
DODGE   <- 0.85                                    # must match position_dodge() below
PAIR_DX <- c(sig = -DODGE / 4, any = DODGE / 4)     # pair centres, offset from the class centre
GAP     <- 0.22 * ymax                              # equal clearance for every pair
ann <- do.call(rbind, lapply(seq_along(ORDER), function(i) {
  cl <- ORDER[i]
  fs <- getf(SCOPE_SIG, cl); fa <- getf(SCOPE_ANY, cl)
  top_sig <- max(getp(SCOPE_SIG, STRATUM, cl)["pct"], getp(SCOPE_SIG, "universe", cl)["pct"])
  top_any <- max(getp(SCOPE_ANY, STRATUM, cl)["pct"], getp(SCOPE_ANY, "universe", cl)["pct"])
  data.frame(
    x   = i + unname(PAIR_DX[c("sig", "any")]),
    y   = unname(c(top_sig, top_any)) + GAP,
    lab = c(sprintf("OR=%.2f\nP=%s", fs["OR"], fmt_p(fs["P"])),
            sprintf("OR=%.2f\nP=%s", fa["OR"], fmt_p(fa["P"]))),
    sig = c(ifelse(fs["P"] < 0.05, "sig", "ns"), ifelse(fa["P"] < 0.05, "sig", "ns")),
    stringsAsFactors = FALSE) }))
Y_TOP <- max(ann$y) + 0.16 * ymax                   # headroom for the two-line labels

fills <- c(SIG, UNI, tint(SIG), tint(UNI)); names(fills) <- G

# subtitle line breaks are the source panel's own; see departure 2 in the header
subtitle_txt <- paste0(
  "Full-strength bars = AS events significant at FDR<0.05 & |dPSI|>=0.10 (as in grouped_pct_overlap_56.pdf); pale bars = ALL AS events detected by rMATS (as in grouped_pct_overlap_56_ANY.pdf).\n",
  "Red = the circRNA group, grey = the whole detected set (universe). The number above each bar counts circRNAs overlapping >=1 event; group sizes (n) are in the legend.\n",
  "Per AS class, a one-sided Fisher's exact test vs universe is annotated for each scope; OR = sample odds ratio.\n",
  "NOTE: the unfiltered scope is saturated, so its ORs sit near 1 (or Inf where no non-matching circRNA remains) and are not evidence of enrichment.")

# ---- build the grouped bar plot (original parameters preserved) ----
p <- ggplot(long, aes(xlab, pct, fill = grp)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = kn), position = position_dodge(width = 0.85), vjust = -0.3, size = 2.0) +
  geom_text(data = ann, aes(x = x, y = y, label = lab), inherit.aes = FALSE, size = 2.4,
            colour = ifelse(ann$sig == "sig", "black", "grey55"), lineheight = 0.9,
            fontface = "bold") +
  scale_fill_manual(values = fills, name = NULL) +
  scale_y_continuous(limits = c(0, Y_TOP), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "AS event class", y = "% of circRNAs sharing \nsplice sites with AS event",
       title = "circRNA x AS overlap by class, both scopes: BSJ-additive-DE (56) vs universe",
       subtitle = subtitle_txt) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title    = element_text(size = 6),                # same size as the subtitle
        plot.subtitle = element_text(size = 6),
        axis.title.x  = element_text(size = axis_render_pt),   # -> 8 pt after the shrink
        axis.text.x   = element_text(size = axis_render_pt),   # -> 8 pt after the shrink
        legend.text   = element_text(size = axis_render_pt)) + # -> 8 pt after the shrink
  # byrow = FALSE fills the 2-row legend COLUMN-wise, so the keys read
  #   top-left dark red / bottom-left dark grey  (HNRNPM-regulated AS)
  #   top-right pale red / bottom-right pale grey (all detected AS)
  # i.e. rows = circRNA set, columns = AS event set. Bar order is set by the
  # factor levels of long$grp and is unaffected by this.
  guides(fill = guide_legend(nrow = 2, byrow = FALSE))

# ---- shared left/right edges for BOTH the legend and the footer ----
# Everything below the panel is squared off against the same two vertical edges:
#   LEFT  = the panel's left edge, which under theme_classic is the y axis itself
#   RIGHT = the right edge of the rightmost bar, which sits INSIDE the panel because the
#           discrete x scale pads 0.6 of a category at each end
# FRAC_R is that right edge as a fraction of the panel width, taken from the BUILT plot so it
# follows the data, the dodge width and the scale expansion instead of being hand-measured.
# Layer 1 is the geom_col layer, so its xmax values are the post-dodge bar edges.
bb     <- ggplot_build(p)
x_rng  <- bb$layout$panel_params[[1]]$x.range          # 0.4 to 6.6 for the six classes
FRAC_R <- (max(bb$data[[1]]$xmax) - x_rng[1]) / diff(x_rng)   # 6.41625 -> 0.9704

# ---- footer: the two group sizes on their own line BELOW the legend ----
# The four legend labels name the circRNA set and the AS event set but carry no counts, so
# the group sizes are stated here instead: LEFT = the red circRNA group (sub_n), RIGHT = the
# grey universe (uni_n). Both are read from the summary table, never hard-coded.
# Why a gtable row and not labs(caption=): ggplot has exactly ONE caption slot with a single
# hjust, so it cannot hold a left-aligned and a right-aligned string on the same line. The
# footer is therefore an extra gtable row added underneath the bottom legend, in the same
# gtable column as the panel, so npc 0 is the y axis and npc FRAC_R is the last bar's edge.
# gtable is a hard dependency of ggplot2, so this adds no new package requirement.
# NOTE: H_CM is unchanged, and the panel is the only height-flexible row, so the panel gives
# up FOOT_ROW of height (about 0.45 cm in the final 16 x 6.272 cm figure) to make room.
FOOT_L   <- sprintf("Total HNRNPM-regulated circRNAs: N = %s", format(sub_n, big.mark = ","))
FOOT_R   <- sprintf("Total number of circRNAs: N = %s",        format(uni_n, big.mark = ","))
FOOT_PT  <- axis_render_pt                     # same rendered size as the legend -> 8 pt final
FOOT_ROW <- grid::unit(FOOT_PT * 1.9, "pt")    # row height = the text plus a little air above
X_L      <- grid::unit(0, "npc")               # the y axis
X_R      <- grid::unit(FRAC_R, "npc")          # right edge of the rightmost bar

# ggplotGrob() freezes the legend's key and label widths into centimetres using the font
# metrics of whatever graphics device happens to be current. Those must be cairo's, or the
# measured label widths differ from the drawn ones by of order 1% and the stretched legend
# overshoots its target edge by a couple of points. So the gtable is built inside a throwaway
# cairo_pdf of exactly the raw output size, matching the device ggsave will use below.
.devfile <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(.devfile, width = W_CM * SC / 2.54, height = H_CM * SC / 2.54)
g <- ggplotGrob(p)
grDevices::dev.off(); unlink(.devfile)

lay <- g$layout
pnl <- lay[lay$name == "panel", ]
gi  <- which(lay$name %in% c("guide-box-bottom", "guide-box"))   # 3.5.x name, then pre-3.5 name
pos <- if (length(gi)) max(lay$b[gi]) else nrow(g) - 1  # else just inside the bottom margin

# ---- stretch the legend to those same two edges ----
# Out of the box ggplot pads the legend with legend.margin and then centres it in its cell
# with a pair of equal null spacers, so it is both narrower than the panel and centred in it.
# Rather than chase the target width through the theme (it is not knowable until the panel
# width is), the legend gtable is lifted out of its guide box and dropped straight into the
# same cell. That cell spans the panel column, so npc inside it IS the panel width. The left
# margin column is zeroed so the first key sits on the y axis, the gap between the legend's
# two columns absorbs whatever slack is needed to reach FRAC_R, and the right margin column
# becomes a filler of the remaining 1 - FRAC_R so the widths sum to exactly 1 npc and no
# justification is required.
#   legend columns: margin | key | label | gap | key | label | margin
# The guard checks that shape (two key columns, the first one second from the left, the last
# label column second from the right) and leaves the legend untouched if it ever differs.
if (length(gi)) {
  LGB <- g$grobs[[gi[1]]]
  li  <- which(LGB$layout$name == "guides")
  if (length(li)) {
    LGD  <- LGB$grobs[[li[1]]]
    LL   <- LGD$layout
    kcol <- sort(unique(LL$l[grepl("^key", LL$name)]))     # key columns, left to right
    lcol <- sort(unique(LL$r[grepl("^label", LL$name)]))   # label columns, left to right
    nW   <- length(LGD$widths)
    if (length(kcol) == 2 && min(kcol) == 2 && max(lcol) == nW - 1) {
      gapc  <- max(kcol) - 1                               # the inter-column gap
      fixed <- sum(LGD$widths[sort(c(kcol, lcol))])         # keys + labels, all absolute
      LGD$widths[1]    <- grid::unit(0, "cm")               # first key flush with the y axis
      LGD$widths[gapc] <- grid::unit(FRAC_R, "npc") - fixed # gap absorbs the slack
      LGD$widths[nW]   <- grid::unit(1 - FRAC_R, "npc")     # filler -> widths sum to 1 npc
      g$grobs[[gi[1]]] <- LGD
    } else {
      warning("legend gtable is not the expected 2-column shape; left as ggplot built it")
    }
  }
}

g <- gtable::gtable_add_rows(g, FOOT_ROW, pos = pos)             # inserts AFTER row `pos`
g <- gtable::gtable_add_grob(
  g, grid::grobTree(
       grid::textGrob(FOOT_L, x = X_L, hjust = 0, gp = grid::gpar(fontsize = FOOT_PT)),
       grid::textGrob(FOOT_R, x = X_R, hjust = 1, gp = grid::gpar(fontsize = FOOT_PT))),
  t = pos + 1, l = min(pnl$l), r = max(pnl$r), clip = "off", name = "group-size-footer")

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
# `g` is a gtable rather than a ggplot, so ggsave cannot read plot.background from the
# theme; bg is given explicitly (theme_classic's plot.background is white anyway).
ggsave(file.path(OUT_DIR, "Fig_S4A_new.pdf"), g,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       bg = "white", device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.4f; raw PDF %.2f x %.2f in; legend/x-axis rendered %.1f pt -> %d pt at %g x %g cm)\n",
            file.path(OUT_DIR, "Fig_S4A_new.pdf"), SC, W_CM * SC / 2.54, H_CM * SC / 2.54,
            axis_render_pt, AXIS_PT, W_CM, H_CM))
cat(sprintf("Footer below the legend: \"%s\"  |  \"%s\"\n", FOOT_L, FOOT_R))

# ============================================================
# FIGURE LEGEND (for Fig_S4A_new.pdf)
# ============================================================
# Co-localisation of HNRNPM-BSJ-DE circRNAs with alternative splicing, broken down by
# alternative-splicing (AS) event class and shown against two definitions of the AS event set.
# For each rMATS event class (SE skipped exon, MXE mutually exclusive exons, A3SS/A5SS
# alternative 3'/5' splice site, RI retained intron, and "Any" = the union across all five
# classes), bars show the percentage of circRNAs whose back-splice splice site coincides with
# the splice site of at least one AS event at the same host gene. Four bars are shown per class.
# The two full-strength bars count only HNRNPM-regulated AS events, called significant at rMATS
# FDR < 0.05 and |dPSI| >= 0.10; the two pale bars count every AS event rMATS detected, with no
# significance filter. Within each pair, red is the group of 56 additive (~ group + PSD)
# BSJ-differentially-expressed circRNAs and grey is the whole detected circRNA universe of
# 5,820; both group sizes are stated in the legend. The number above each bar is the count of
# circRNAs in that group overlapping at least one such event. Above each pair, the sample odds
# ratio (OR) and the one-sided Fisher's exact test P value (alternative, group more enriched
# than universe) are annotated in bold, printed in grey where P >= 0.05, and each label is
# placed the same fixed distance above the taller bar of its own pair. The horizontal axis is
# the AS event class and the vertical axis is the percentage of circRNAs overlapping an AS
# event. BSJ-DE circRNAs are enriched for co-localising with HNRNPM-differential AS events
# overall and specifically at skipped-exon (SE) and mutually-exclusive-exon (MXE) events. The
# unfiltered comparison is included as a control and not as evidence: because 95.6% of the
# universe already sits at the splice site of some rMATS-tested event, that comparison is
# saturated, its odds ratios sit near 1, and it reports OR = Inf wherever a group contains no
# non-overlapping circRNA (a zero denominator, not an infinitely strong effect — the
# accompanying P value is non-significant in those cells). P values are per class and per
# scope and are not corrected for multiple testing.
