# Fig S5C — proliferation-only circRNA fold-change ceiling, gHNRNPM_1, Day 4 S-phase.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w3_circ_FC_ceiling/doubling_vs_genotype.R
# The source draws a 2-page PDF; page 2 ("Supplementary", Day 4) is
#   panel(gHNRNPM_1,Day 4) | panel(gHNRNPM_2,Day 4). This reproduces ONLY the gHNRNPM_1 Day 4
#   panel. Day-4 counterpart of Fig_5E.R. Rather than re-running the growth-window selection,
#   it reads the scalars the source already wrote to manuscript_inputs.csv and rebuilds the
#   panel with the identical model curve, ceilings and annotations. (The growth-curve ceiling
#   is timepoint-independent, so only the EdU ceiling ratio differs from Day 1 -> R_edu_Day4.)
#
# Strategy (same as the Fig_2/Fig_3/Fig_4/Fig_5 series): keep ALL original plot parameters,
#   render at the original size, and let ggsave(scale=) squash the figure to 7.75 x 6.83 cm
#   (uniform shrink preserves proportions). Pre-inflate ONLY the axis title/tick fonts + legend
#   so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger than
#   nominal; resize it to 7.75 cm wide on the slide. See Fig_2E_exactsize.R.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
W3_DIR  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w3_circ_FC_ceiling"
INPUTS  <- file.path(W3_DIR, "manuscript_inputs.csv")
OUT_DIR <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/Figures/Supplementary"

GUIDE <- "gHNRNPM_1"; DAY <- "Day 4"

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 6                     # one panel of the source page (12 in wide, 2 panels) is ~6 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

ln2 <- log(2)

# ---- read the scalars the source computed (single source of truth) ----
minp <- read.csv(INPUTS, stringsAsFactors = FALSE, check.names = FALSE)
getv <- function(metric, guide) as.numeric(minp$value[minp$metric == metric & minp$guide == guide])
gets <- function(metric, guide) as.character(minp$value[minp$metric == metric & minp$guide == guide])

Redu    <- getv("R_edu_Day4", GUIDE)   # EdU-derived ceiling ratio at Day 4
Rgro    <- getv("R_growth",   GUIDE)   # growth-curve doubling-time ceiling ratio
Tctl    <- getv("Td_ctrl_h",  GUIDE)   # control doubling time (hours)
win_lab <- gets("win_label",  GUIDE)   # chosen exponential window, e.g. "D0-D3"
FC_min  <- getv("FC_min", "all")       # smallest independent-circRNA fold change
n_circ  <- as.integer(getv("n_circ", "all"))

# ---- proliferation-only fold-change model (identical to source) ----
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

half <- 10 ^ seq(log10(3), log10(1000), length.out = 400)
df <- rbind(
  data.frame(half = half, fc = prolif_fc(half, Tctl, Redu), src = "EdU S-phase ratio"),
  data.frame(half = half, fc = prolif_fc(half, Tctl, Rgro), src = "Growth-curve doubling time"))

y.top <- max(FC_min, Redu, Rgro) * 1.06

# growth-curve ceiling label: place BELOW its dashed line when the growth ceiling is
# the lower of the two ratios, ABOVE when it is the higher one (keeps it clear of the EdU label)
gro_vjust <- if (Rgro < Redu) 1.4 else -0.6

# ---- build the panel (original parameters preserved) ----
pg <- ggplot(df, aes(half, fc, colour = src)) +
  geom_hline(yintercept = Redu, linetype = "dashed", colour = col_edu,    linewidth = 0.6) +
  geom_hline(yintercept = Rgro, linetype = "dashed", colour = col_growth, linewidth = 0.6) +
  geom_hline(yintercept = FC_min, colour = col_circ, linetype = "dotted", linewidth = 0.7) +
  geom_line(linewidth = 1) +
  annotate("text", x = 3, y = FC_min, hjust = 0, vjust = -0.5, size = 4.5, colour = col_circ,
           label = paste0("smallest independent circRNA FC = ", round(FC_min, 2),
                          "  (n = ", n_circ, ")")) +
  annotate("text", x = 3, y = Redu, hjust = 0, vjust = -0.6, size = 4.5, colour = col_edu,
           label = paste0("EdU ceiling ratio = ", round(Redu, 2))) +
  annotate("text", x = 3, y = Rgro, hjust = 0, vjust = gro_vjust, size = 4.5, colour = col_growth,
           label = paste0("growth-curve ceiling ratio = ", round(Rgro, 2))) +
  scale_colour_manual(values = c("EdU S-phase ratio" = col_edu,
                                 "Growth-curve doubling time" = col_growth),
                      name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +   # stack the two legend keys vertically
  scale_x_log10(breaks = c(3, 10, 30, 100, 300, 1000),
                labels = c("3", "10", "30", "100", "300", "1000")) +
  coord_cartesian(ylim = c(1.0, y.top)) +
  labs(title = paste0(GUIDE, "  (", DAY, " S-phase)"),
       subtitle = paste0("growth-curve exponential window ", win_lab),
       x = "circRNA half-life (hours)",
       y = "circRNA fold increase (KD/Ctrl) attributable \nto reduced proliferation in HNRNPM-KD") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title  = element_text(face = "bold"),     # original proportion
        axis.title  = element_text(size = axis_render_pt),   # -> 8 pt after shrink
        axis.text   = element_text(size = 0.5 * axis_render_pt),   # -> 4 pt after shrink (50% of previous ticks)
        legend.text = element_text(size = 0.5 * axis_render_pt),   # -> 4 pt after shrink (50% of previous)
        legend.key.spacing.y = grid::unit(0, "pt"),   # minimal gap between the two stacked legend keys
        legend.box.spacing   = grid::unit(0, "pt"),   # minimal gap between legend and x-axis label
        legend.margin        = margin(0, 0, 0, 0))

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_S5C.pdf"), pg,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_S5C.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_S5C.pdf)
# ============================================================
# Proliferation-only ceiling on the circRNA fold increase that reduced division could
# produce in gHNRNPM_1 knockdown, using the Day 4 EdU data. Under a steady-state
# single-compartment model of a stable RNA, a slower-dividing cell dilutes the RNA less,
# so the proliferation-only fold change (knockdown / control) is
# (k_deg + mu_ctrl) / (k_deg + mu_KD), with k_deg = ln2/half-life and mu = ln2/doubling
# time; as half-life increases this approaches its ceiling, the knockdown-to-control
# doubling-time ratio R. The two coloured curves show the modelled fold change versus
# circRNA half-life for R estimated two ways: from the EdU S-phase ratio (blue) and from
# the growth-curve doubling-time ratio over the objectively chosen exponential window
# (orange); the matching dashed lines mark each ceiling. The dotted black line is the
# smallest fold change among the 36 HNRNPM-independent circRNAs (Supplementary Table
# S5B, logFC_BSJ). Because even that smallest circRNA fold change exceeds both ceilings,
# reduced proliferation alone cannot account for the circRNA increases. The horizontal
# axis is circRNA half-life in hours (log scale) and the vertical axis is the
# proliferation-attributable fold increase.
