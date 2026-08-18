# Fig S5A — EdU % S-phase bar plot, gHNRNPM_1 vs gCtrl, Day 4.
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w1_EdU/1_EdU_statistics_twotests.R
# The source draws four comparisons (gHNRNPM_1/2 x Day 1/4) as four PDF pages; this
#   reproduces ONLY the gHNRNPM_1 Day 4 panel, recomputing the Student's and Welch's
#   two-sample t-tests exactly as the source and drawing the same annotated bar plot.
#   Sister of Fig_5C.R (gHNRNPM_1 Day 1); Day 4 uses cols D/E.
#
# Strategy (same as the Fig_2/Fig_3/Fig_4/Fig_5 series): keep ALL original plot
#   parameters (bar width/spacing, jitter points, significance bracket, p-value
#   annotation, title/subtitle/caption sizes), render at the original size, and let
#   ggsave(scale=) squash the figure to 7.75 x 6.83 cm (uniform shrink preserves
#   proportions). Pre-inflate ONLY the axis title/tick fonts so they land at 8 pt.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x larger
#   than nominal; resize it to 7.75 cm wide on the slide. See Fig_2E_exactsize.R.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(readxl); library(ggplot2)
})


# --- Path definition ---
EDU_DIR    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w1_EdU"
INPUT_XLSX <- file.path(EDU_DIR, "Raw_percentage.xlsx")
SHEET_NAME <- "% S-phase"
OUT_DIR    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/Figures/Supplementary"

ALTERNATIVE <- "two.sided"
BAR_WIDTH   <- 0.45
BAR_SPACING <- 0.6
BASE_FONT   <- 12

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 85 / 25.4             # source page width = FIG_WIDTH_MM/25.4 = 85 mm
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink

# ---- Read the sheet; gHNRNPM_1 Day 4 is cols D/E (4/5), rows 5-7 ----
raw <- as.data.frame(read_excel(INPUT_XLSX, sheet = SHEET_NAME,
                                col_names = FALSE, .name_repair = "minimal"))
getvals <- function(rows, col) { v <- suppressWarnings(as.numeric(raw[rows, col])); v[!is.na(v)] }
guide <- "gHNRNPM_1"; day <- "Day 4"
ctrl  <- getvals(5:7, 4)
treat <- getvals(5:7, 5)

# ---- helpers (identical to source) ----
sem   <- function(x) sd(x) / sqrt(length(x))
stars <- function(p) if (is.na(p)) "n/a" else if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < 5e-2) "*" else "ns"
p_expr <- function(p) {
  if (is.na(p)) return("italic(p)*' = NA'")
  exponent <- floor(log10(p)); mantissa <- round(p / 10^exponent, 2)
  if (mantissa >= 10) { mantissa <- mantissa / 10; exponent <- exponent + 1 }
  paste0("italic(p)=='", formatC(mantissa, format = "f", digits = 2), "'%*%10^", exponent)
}

# ---- both t-tests (Student's shown above bars, Welch's in caption) ----
tt_s <- t.test(ctrl, treat, alternative = ALTERNATIVE, var.equal = TRUE)
tt_w <- t.test(ctrl, treat, alternative = ALTERNATIVE, var.equal = FALSE)
p_s  <- tt_s$p.value; p_w <- tt_w$p.value

# ---- tidy frames (bars at explicit x-positions) ----
grp_levels <- c("gCtrl", guide)
xpos       <- c(1, 1 + BAR_SPACING)
raw_df <- data.frame(group = factor(rep(grp_levels, c(length(ctrl), length(treat))), levels = grp_levels),
                     x = rep(xpos, c(length(ctrl), length(treat))), value = c(ctrl, treat))
summ_df <- data.frame(group = factor(grp_levels, levels = grp_levels), x = xpos,
                      mean = c(mean(ctrl), mean(treat)), sem = c(sem(ctrl), sem(treat)))

ymax   <- max(raw_df$value, summ_df$mean + summ_df$sem)
bar_y  <- ymax * 1.10; tick <- ymax * 0.02; text_y <- bar_y + ymax * 0.06
cols <- c("#9E9E9E", "#D55E00"); names(cols) <- grp_levels
ann_size <- BASE_FONT / .pt                  # original annotation text size (pt -> mm)
welch_cap <- parse(text = paste0("\"Welch's two-sample t-test: \"*", p_expr(p_w)))

# ---- build the bar plot (original parameters preserved) ----
pg <- ggplot(summ_df, aes(x = x, y = mean, fill = group)) +
  geom_col(width = BAR_WIDTH, colour = "black") +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = BAR_WIDTH * 0.4) +
  geom_jitter(data = raw_df, aes(x = x, y = value), width = BAR_WIDTH * 0.25, height = 0,
              size = 2, shape = 21, fill = "white", colour = "black", inherit.aes = FALSE) +
  annotate("segment", x = xpos[1], xend = xpos[2], y = bar_y, yend = bar_y) +
  annotate("segment", x = xpos[1], xend = xpos[1], y = bar_y, yend = bar_y - tick) +
  annotate("segment", x = xpos[2], xend = xpos[2], y = bar_y, yend = bar_y - tick) +
  # significance asterisk only (Student's p); numeric p omitted since the star conveys it
  annotate("text", x = mean(xpos), y = text_y,
           label = paste0("'", stars(p_s), "'"), parse = TRUE,
           size = ann_size, fontface = "bold") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(breaks = xpos, labels = grp_levels, expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = paste0(guide, " vs gCtrl - ", day),
       subtitle = "Student's two-sample t-test",
       x = NULL, y = "% cells in S-phase", caption = welch_cap) +
  theme_classic(base_size = BASE_FONT) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = BASE_FONT),  # original proportion
        plot.subtitle = element_text(hjust = 0.5, size = BASE_FONT),                 # original proportion
        plot.caption  = element_text(hjust = 0.5, size = BASE_FONT - 2, colour = "black"),  # original proportion
        axis.text     = element_text(size = axis_render_pt),                         # -> 8 pt after shrink
        axis.title    = element_text(size = axis_render_pt),                         # -> 8 pt after shrink
        axis.text.x   = element_text(face = "bold", size = axis_render_pt))          # -> 8 pt after shrink (bold)

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_S5A.pdf"), pg,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (scale=%.3f; axis rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_S5A.pdf"), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_S5A.pdf)
# ============================================================
# EdU incorporation (percentage of cells in S-phase) four days after CRISPRi knockdown
# of HNRNPM with guide gHNRNPM_1, compared to a non-targeting control (gCtrl), in
# SH-SY5Y (dCas9-KRAB) cells. Bars show the group mean, error bars the standard error
# of the mean, and open circles the individual replicate wells (gCtrl grey, gHNRNPM_1
# orange). The horizontal axis is the guide group and the vertical axis is the
# percentage of cells in the S-phase gate. The bracket above the bars marks the
# significance of the Student's (equal-variance) two-sample two-sided t-test with
# asterisks (* P < 0.05, ** P < 0.01, *** P < 0.001, ns not significant); here the
# Student's two-sample two-sided t-test gives P = 1.92 x 10^-2 (*). The caption reports
# the Welch's (unequal-variance) two-sample t-test P value on the same data.
