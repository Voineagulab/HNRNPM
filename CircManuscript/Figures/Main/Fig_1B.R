# Fig 1B — HnrnpM eCLIP enrichment at circRNA flanking introns (box plot only).
# Adapted from: /mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/SCRIPTS/09_HnrnpM_CLIP_enrichment/plot_enrichment_PSD58.R
# Retains only PANEL A of that figure — the per-region log2(IP/Input) box plot
#   (the upper panel; the lower violin+box panel B is dropped).
# DE = additive BSJ-DE circRNAs (~ group + PSD, adj.P<0.05, publication primary).
# Saved as PDF sized for the manuscript slide (W_CM x H_CM), fonts 8 pt, via
#   cairo_pdf, to Figures/Main/Fig_1B.pdf.

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
})


# --- Path definition ---
ROOT    <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV"
RES09   <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/09_HnrnpM_CLIP_enrichment")
SIG     <- file.path(RES09, "clip_signal_per_region.tsv")
BSJ_ADD <- file.path(ROOT, "ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/limma_BSJ_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
OUT_DIR <- file.path(ROOT, "JW/CircManuscript/Figures/Main")

# ---- Target physical size on the slide ----
W_CM   <- 5.77; H_CM <- 6.83      # width x height as it should appear in PowerPoint
MIN_PT <- 8                        # font size (points)

# ---- Load data ----
sig <- read.table(SIG,     sep = "\t", header = TRUE, check.names = FALSE)
bsj <- read.table(BSJ_ADD, sep = "\t", header = TRUE, check.names = FALSE)

# Additive BSJ-DE circRNAs (publication primary)
de_circs <- bsj$circRNA[bsj$adj.P.Val < 0.05]

# ---- Per-region log2(IP/Input) box plot (panel A of the source figure) ----
sig$de_lab <- ifelse(sig$circRNA %in% de_circs, "DE", "non-DE")
sig$panel_label <- with(sig, ifelse(
  region_type == "flank5", paste0(de_lab, " 5' flank"),
  ifelse(region_type == "flank3", paste0(de_lab, " 3' flank"),
  ifelse(region_type == "within_gene_bg", "Within-gene bg",
  ifelse(region_type == "non_host_bg",   "Non-host-gene bg", region_type)))))
sig$panel_label <- factor(sig$panel_label,
  levels = c("DE 5' flank","DE 3' flank","non-DE 5' flank","non-DE 3' flank",
             "Within-gene bg","Non-host-gene bg"))
palette <- c("DE 5' flank" = "#d73027", "DE 3' flank" = "#d73027",
             "non-DE 5' flank" = "#4C78A8", "non-DE 3' flank" = "#4C78A8",
             "Within-gene bg" = "grey50", "Non-host-gene bg" = "grey75")

pb <- ggplot(sig, aes(x = panel_label, y = log2_IP_over_Input, fill = panel_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, colour = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = palette, guide = "none") +
  labs(title = "log2(IP/Input) per region",
       x = NULL, y = "log2(IP / Input)") +
  theme_bw(base_size = MIN_PT) +
  theme(
    text        = element_text(size = MIN_PT),
    axis.text   = element_text(size = MIN_PT),
    axis.title  = element_text(size = MIN_PT),
    plot.title  = element_text(size = MIN_PT, hjust = 0.5, margin = margin(b = 1)),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(OUT_DIR, "Fig_1B.pdf"), pb,
       width = W_CM, height = H_CM, units = "cm",
       device = cairo_pdf)
cat("Wrote:", file.path(OUT_DIR, "Fig_1B.pdf"), "\n")


# ---- Figure legend ----
# (B) Box plot of HNRNPM eCLIP enrichment across six genomic region categories. eCLIP signal is expressed as log2 of immunoprecipitation coverage over input coverage (IP/Input), computed per region from HNRNPM eCLIP data published by Ho et al (2021) lifted to the hg38 assembly. The categories are the 5 prime and 3 prime introns flanking the 56 DE circRNAs highlighted in (A), the 5 prime and 3 prime introns flanking non-DE circRNAs, a within DE circRNA host gene intronic background, and a background drawn from introns of expressed genes that do not host DE circRNAs. Boxes show the median and interquartile range with whiskers, and the horizontal dashed line marks equal immunoprecipitation and input signal at a log2 value of zero. Higher values indicate stronger HNRNPM binding.