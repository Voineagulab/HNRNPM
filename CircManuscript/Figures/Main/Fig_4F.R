# Fig 4F — rMATS hit counts, PSD5+8 (full) run only.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/07_Splicing_rMATS/summarize_rmats_PSD58.R
# Reproduces ONLY the PSD5+8 (full) hit-count barplot from `rmats_summary_combined_PSD58.pdf`
#   (the top panel's PSD5+8 facet). The PSD5-only facet and the overlap barplot are dropped.
#   Reads the saved full-run summary table directly.
#
# Strategy (same as the Fig_2/Fig_3/Fig_4 series): keep ALL original plot parameters,
#   render at the original size, and let ggsave(scale=) squash the figure to
#   5.46 x 6.83 cm (uniform shrink preserves proportions). Pre-inflate ONLY the axis
#   title/tick fonts + legend so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger
#   than nominal; resize it to 5.46 cm wide on the slide. See Fig_2E_exactsize.R.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES07    <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/07_Splicing_rMATS")
FULL_SUM <- file.path(RES07, "rmats_PSD58_full_summary.tsv")
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

FDR_THR <- 0.05; DPSI_THR <- 0.10   # rMATS significance thresholds (for the title)

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 5.46; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 5.5                   # one facet of the source top panel (11 in wide, 2 facets) is ~5.5 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Load full-run summary; build hit-count long table (PSD5+8 only) ----
full <- read.table(FULL_SUM, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
hc <- rbind(
  data.frame(event_type = full$event_type, direction = "Inclusion UP",   n = full$n_up),
  data.frame(event_type = full$event_type, direction = "Inclusion DOWN", n = full$n_down))

p_hc <- ggplot(hc, aes(x = event_type, y = n, fill = direction)) +
  geom_col(position = "dodge", width = 0.8) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Inclusion UP" = "#d73027", "Inclusion DOWN" = "#2166ac"),
                    labels = c("Inclusion UP" = "UP", "Inclusion DOWN" = "DOWN")) +
  labs(x = NULL, y = "# sig events", fill = "") +
  theme_bw(base_size = 11) +                          # original base size / proportions
  theme(
    axis.title      = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    axis.text       = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.title    = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.text     = element_text(size = axis_render_pt),  # -> 8 pt after the ggsave shrink
    legend.position = "bottom"                              # legend moved from right to bottom
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_4F.pdf"), p_hc,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_4F.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_4F.pdf)
# ============================================================
# Differential alternative splicing between control and HNRNPM knockdown in the
# PSD5+8 contrast (3 knockdown vs 4 control), from rMATS-turbo. Bars show the number
# of significant events per event type (SE skipped exon, MXE mutually exclusive
# exons, A3SS/A5SS alternative 3'/5' splice site, RI retained intron), called
# significant at FDR < 0.05 and |dPSI| >= 0.10, split by the direction of the
# inclusion change in HNRNPM knockdown (inclusion up in red, inclusion down in blue).
# The count for each bar is printed above it. The horizontal axis is the event type
# and the vertical axis is the number of significant events.
