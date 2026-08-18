# Fig 4C — GO Cellular Component enrichment dot plot, HNRNPM-KD UP-regulated genes.
# Adapted from the `go_dotplot()` block of
#   /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/Manuscript_R/SCRIPTS/render_new_figures.R
#   (which rebuilds panel_GO_UP_CC.png), using the FinalAnalysis enrichment table
#   HnrnpM_UP_GO_CC.tsv directly. Top 20 enriched GO CC terms as a dot plot
#   (x = GeneRatio, y = term, size = Count, colour = p.adjust).
#
# Strategy (same as the Fig_2/Fig_3/Fig_4 series): keep ALL original plot parameters
#   (dot size range, colour gradient, title/subtitle sizes, wrapped term labels),
#   render at the original size, and let ggsave(scale=) squash the figure to
#   7.75 x 6.83 cm (uniform shrink preserves proportions). Pre-inflate ONLY the axis
#   title/tick fonts + legend so they land at 8 pt in the output.
# NOTE: ggsave(scale=) multiplies the OUTPUT dimensions, so the raw PDF is SC x
#   larger than nominal; resize it to 7.75 cm wide on the slide (that lands the axis
#   at 8 pt). See Fig_2E_exactsize.R for the born-at-size variant.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT     <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
ENR_DIR  <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/02_GeneTx_DE/Enrichment")
GO_TSV   <- file.path(ENR_DIR, "HnrnpM_UP_GO_CC.tsv")
OUT_DIR  <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size + uniform-shrink bookkeeping ----
W_CM   <- 7.75; H_CM <- 6.83      # nominal final PDF size on the slide
ORIG_W_IN <- 9                     # the source rendered this dot plot at width = 9 in
SC     <- (ORIG_W_IN * 2.54) / W_CM    # ggsave scale: draw at original size, shrink to target
AXIS_PT <- 8                       # desired FINAL axis + legend font size
axis_render_pt <- AXIS_PT * SC     # pre-inflate so it becomes 8 pt after the 1/SC shrink
GO_TOP_N <- 20

# ---- Load enrichment table, take top 20 by p.adjust ----
d0 <- read.table(GO_TSV, sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                 stringsAsFactors = FALSE)
n_total <- nrow(d0)                                        # all significant terms (BH p.adjust<0.05)
d0$gr <- vapply(strsplit(d0$GeneRatio, "/"),
                function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
d <- d0[order(d0$p.adjust), ]
d <- d[seq_len(min(GO_TOP_N, n_total)), ]
wrap_lab <- function(x, w = 40) vapply(x, function(s) paste(strwrap(s, width = w), collapse = "\n"), character(1))
d$lab <- wrap_lab(d$Description)
d$lab <- factor(d$lab, levels = d$lab[order(d$gr)])        # highest GeneRatio at top

# ---- Dot plot (original parameters preserved) ----
g <- ggplot(d, aes(x = gr, y = lab, size = Count, colour = p.adjust)) +
  geom_point() +
  scale_colour_gradient(low = "#d73027", high = "#4575b4", name = "p.adjust") +
  scale_size_continuous(range = c(2, 7), name = "Count") +
  labs(title = "HNRNPM-KD UP-regulated genes: GO Cellular Component",
       subtitle = sprintf("top %d of %d enriched terms (BH p.adjust < 0.05)", nrow(d), n_total),
       x = "GeneRatio", y = NULL) +
  theme_bw(base_size = 16.5) +                        # original base size / proportions
  theme(
    axis.text.y   = element_text(size = 0.6 * axis_render_pt, lineheight = 0.7),  # -> 4.8 pt (60%); tighter line spacing for wrapped terms
    axis.text.x   = element_text(size = axis_render_pt),                     # -> 8 pt after shrink
    axis.title    = element_text(size = axis_render_pt),                     # -> 8 pt after shrink
    legend.title  = element_text(size = axis_render_pt),                     # -> 8 pt after shrink
    legend.text   = element_text(size = axis_render_pt),                     # -> 8 pt after shrink
    plot.title    = element_text(face = "bold", size = 18),                  # original proportion
    plot.subtitle = element_text(size = 13.5, colour = "grey40"),           # original proportion
    panel.grid.minor = element_blank()
  )

# scale = SC renders at the original size, then the file is squashed to W_CM x H_CM.
ggsave(file.path(OUT_DIR, "Fig_4C.pdf"), g,
       width = W_CM, height = H_CM, units = "cm", scale = SC,
       device = cairo_pdf)
cat(sprintf("Wrote: %s (%d terms; scale=%.3f; axis/legend rendered %.1f pt -> %d pt final)\n",
            file.path(OUT_DIR, "Fig_4C.pdf"), nrow(d), SC, axis_render_pt, AXIS_PT))

# ============================================================
# FIGURE LEGEND (for Fig_4C.pdf)
# ============================================================
# Gene Ontology Cellular Component (GO CC) over-representation analysis of the genes
# UP-regulated in HNRNPM knockdown (additive ~ group + PSD gene-level DE). Each dot
# is one enriched GO CC term; the top 20 terms by adjusted P value are shown, ordered
# with the highest GeneRatio at the top. The horizontal axis is the GeneRatio (the
# fraction of the queried UP gene set annotated to that term), dot size encodes the
# number of UP genes in the term (Count), and dot colour encodes the BH-adjusted
# enrichment P value (p.adjust). Only terms significant at BH p.adjust < 0.05 are
# included; the subtitle reports how many of the total significant terms are shown.
