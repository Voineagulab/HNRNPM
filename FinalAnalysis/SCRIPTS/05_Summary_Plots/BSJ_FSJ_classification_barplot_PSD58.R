suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })

# Summary bar plot for the BSJ x FSJ mechanistic classification of additive
# BSJ-DE circRNAs (publication primary). Adds the primary-model version as
# a small comparator.

ROOT <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ"
PRIMARY_TBL <- file.path(ROOT, "BSJ_FSJ_classification.tsv")
ADD_TBL     <- file.path(ROOT, "BSJ_FSJ_classification_additive.tsv")
OUT_PDF <- file.path("/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/05_Summary_Plots",
                     "BSJ_FSJ_classification_barplot_PSD58.pdf")

pri <- read.table(PRIMARY_TBL, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
add <- read.table(ADD_TBL,     sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Restrict to BSJ-DE circRNAs in each (note: primary table uses 'sig_BSJ' bool from primary
# DE; additive table uses 'sig_BSJ' bool from additive DE)
to_bool <- function(x) x == "True" | x == TRUE
pri$sig_BSJ <- to_bool(pri$sig_BSJ)
add$sig_BSJ <- to_bool(add$sig_BSJ); add$fsj_filter_fail <- to_bool(add$fsj_filter_fail)
# Primary table (built by limma_voom_FSJ_PSD58.R) already encodes linear-collapse
# via the "extreme_linear_collapse" class derived from the FSJ filter pass status.
# The additive table separately carries fsj_filter_fail. Not needed downstream
# because both share the `class` column.

pri_sig <- pri[pri$sig_BSJ, ]
add_sig <- add[add$sig_BSJ, ]

# Re-tabulate counts robustly (use the 'class' column which already encodes the 4-way)
class_levels <- c("independent_backsplicing","co_regulated","opposite_regulation","extreme_linear_collapse")
counts <- data.frame(
  class = factor(class_levels, levels = class_levels),
  primary  = sapply(class_levels, function(k) sum(pri_sig$class == k)),
  additive = sapply(class_levels, function(k) sum(add_sig$class == k))
)
counts$primary_total  <- sum(counts$primary)
counts$additive_total <- sum(counts$additive)
counts$primary_pct  <- round(100 * counts$primary  / counts$primary_total,  1)
counts$additive_pct <- round(100 * counts$additive / counts$additive_total, 1)

long <- rbind(
  data.frame(model = sprintf("Additive (~group+PSD, n=%d)", counts$additive_total[1]),
             class = counts$class, n = counts$additive, pct = counts$additive_pct),
  data.frame(model = sprintf("Primary (~group, n=%d)", counts$primary_total[1]),
             class = counts$class, n = counts$primary,  pct = counts$primary_pct))
long$model <- factor(long$model, levels = c(sprintf("Additive (~group+PSD, n=%d)", counts$additive_total[1]),
                                              sprintf("Primary (~group, n=%d)",      counts$primary_total[1])))

pal <- c("independent_backsplicing"="#1f77b4","co_regulated"="#2ca02c",
         "opposite_regulation"="#d62728","extreme_linear_collapse"="#9467bd")
p <- ggplot(long, aes(x = class, y = pct, fill = class)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)), vjust = -0.4, size = 3) +
  facet_wrap(~ model, nrow = 1) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "BSJ × FSJ mechanistic class of HnrnpM-BSJ-DE circRNAs",
       subtitle = "Publication primary = additive model. Independent back-splicing dominates under both designs.",
       x = NULL, y = "% of BSJ-DE circRNAs") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey92", colour = "grey60"))
ggsave(OUT_PDF, p, width = 11, height = 6)
cat("Wrote:", OUT_PDF, "\n")
