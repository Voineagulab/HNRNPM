suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# FinalAnalysis 07 — summarize the two rMATS runs at FDR<0.05 & |dPSI|>=0.10:
#   PSD5+8  (full 7-sample cohort, 3 HnrnpM + 4 NEG4) — primary splicing result
#   PSD5-only (4 samples, 2v2)                       — early-KD sensitivity
# Writes per-event-type hit counts, the combined sig event table for each run,
# and a coordinate-based overlap table between the two runs.

ROOT       <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis"
RAW_FULL   <- file.path(ROOT, "RESULTS/07_Splicing_rMATS/rmats_raw")
RAW_PSD5   <- file.path(ROOT, "RESULTS/07_Splicing_rMATS/rmats_PSD5_only_HnrnpM")
OUT_DIR    <- file.path(ROOT, "RESULTS/07_Splicing_rMATS")
EVENT_TYPES <- c("SE","RI","A3SS","A5SS","MXE")
FDR_THR  <- 0.05; DPSI_THR <- 0.10

COORD_COLS <- list(
  SE   = c("exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"),
  RI   = c("riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"),
  A3SS = c("longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE"),
  A5SS = c("longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE"),
  MXE  = c("1stExonStart_0base","1stExonEnd","2ndExonStart_0base","2ndExonEnd",
           "upstreamES","upstreamEE","downstreamES","downstreamEE"))

read_event <- function(et, dir) {
  f <- file.path(dir, paste0(et, ".MATS.JCEC.txt"))
  if (!file.exists(f)) { cat("Missing:", f, "\n"); return(NULL) }
  df <- read.table(f, sep = "\t", header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "")
  df$event_type <- et
  df
}

summarise_one <- function(df, et, label) {
  if (is.null(df)) return(NULL)
  sig <- df$FDR < FDR_THR & !is.na(df$FDR) &
         abs(df$IncLevelDifference) >= DPSI_THR & !is.na(df$IncLevelDifference)
  data.frame(run = label, event_type = et, n_tested = nrow(df),
             n_sig = sum(sig),
             n_up = sum(sig & df$IncLevelDifference > 0),
             n_down = sum(sig & df$IncLevelDifference < 0))
}

event_key <- function(df, et) {
  cols <- COORD_COLS[[et]]
  paste(df$GeneID, df$chr, df$strand,
        apply(df[, cols, drop = FALSE], 1, paste, collapse = "_"), sep = "|")
}

keep_cols <- c("ID","GeneID","geneSymbol","chr","strand",
               "PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference","event_type")

# Per-run: per-event-type sig table + combined sig table
collect <- function(raw_dir, label) {
  events <- lapply(EVENT_TYPES, read_event, dir = raw_dir); names(events) <- EVENT_TYPES
  summ <- do.call(rbind, Map(summarise_one, events, EVENT_TYPES, label))
  sig  <- do.call(rbind, lapply(events, function(d) {
    if (is.null(d)) return(NULL)
    keep <- d$FDR < FDR_THR & !is.na(d$FDR) &
            abs(d$IncLevelDifference) >= DPSI_THR & !is.na(d$IncLevelDifference)
    if (!any(keep)) return(NULL)
    cp <- intersect(keep_cols, colnames(d))
    d[keep, cp, drop = FALSE]
  }))
  list(events = events, summ = summ, sig = sig)
}

full_run <- collect(RAW_FULL, "PSD5+8 full (3HM vs 4NEG)")
psd5_run <- collect(RAW_PSD5, "PSD5-only (2HM vs 2NEG)")

# Save per-event summaries + combined sig tables
write.table(full_run$summ, file.path(OUT_DIR, "rmats_PSD58_full_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(full_run$sig,  file.path(OUT_DIR, "rmats_PSD58_full_sig.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(psd5_run$summ, file.path(OUT_DIR, "rmats_PSD5_only_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(psd5_run$sig,  file.path(OUT_DIR, "rmats_PSD5_only_sig.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== PSD5+8 full (publication primary splicing result) ===\n")
print(full_run$summ, row.names = FALSE)
cat("\n=== PSD5-only (early-KD sensitivity) ===\n")
print(psd5_run$summ, row.names = FALSE)

# Coordinate-based overlap PSD5+8 vs PSD5-only
overlap_rows <- list()
for (et in EVENT_TYPES) {
  full_d <- full_run$events[[et]]; psd5_d <- psd5_run$events[[et]]
  if (is.null(full_d) || is.null(psd5_d)) next
  full_d$key <- event_key(full_d, et); psd5_d$key <- event_key(psd5_d, et)
  full_sig_mask <- !is.na(full_d$FDR) & full_d$FDR < FDR_THR &
                   !is.na(full_d$IncLevelDifference) & abs(full_d$IncLevelDifference) >= DPSI_THR
  psd5_sig_mask <- !is.na(psd5_d$FDR) & psd5_d$FDR < FDR_THR &
                   !is.na(psd5_d$IncLevelDifference) & abs(psd5_d$IncLevelDifference) >= DPSI_THR
  full_sig_keys <- full_d$key[full_sig_mask]; psd5_sig_keys <- psd5_d$key[psd5_sig_mask]
  shared <- length(intersect(full_sig_keys, psd5_sig_keys))
  in_psd5_tested  <- sum(full_sig_keys %in% psd5_d$key)
  if (in_psd5_tested > 0) {
    pair <- merge(
      data.frame(key = full_sig_keys, dPSI_full = full_d$IncLevelDifference[full_sig_mask], stringsAsFactors = FALSE),
      data.frame(key = psd5_d$key,    dPSI_psd5 = psd5_d$IncLevelDifference, stringsAsFactors = FALSE),
      by = "key")
    same_dir <- sum(sign(pair$dPSI_full) == sign(pair$dPSI_psd5), na.rm = TRUE)
  } else same_dir <- 0
  overlap_rows[[et]] <- data.frame(
    event_type = et,
    n_full_sig = length(full_sig_keys), n_PSD5_sig = length(psd5_sig_keys),
    n_shared_sig = shared,
    pct_full_in_PSD5_sig = ifelse(length(full_sig_keys)>0,
                                  round(100 * shared / length(full_sig_keys), 1), NA),
    n_full_in_PSD5_tested = in_psd5_tested,
    n_same_dir_in_PSD5 = same_dir,
    pct_same_dir_among_tested = ifelse(in_psd5_tested > 0,
                                       round(100 * same_dir / in_psd5_tested, 1), NA))
}
overlap <- do.call(rbind, overlap_rows)
write.table(overlap, file.path(OUT_DIR, "rmats_PSD5_only_vs_full_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Overlap: PSD5-only vs PSD5+8 full ===\n"); print(overlap, row.names = FALSE)

# Hit-count + overlap bar plots
hc <- rbind(
  data.frame(run = "PSD5+8 (full)",    event_type = full_run$summ$event_type,
             direction = "Inclusion UP", n = full_run$summ$n_up),
  data.frame(run = "PSD5+8 (full)",    event_type = full_run$summ$event_type,
             direction = "Inclusion DOWN", n = full_run$summ$n_down),
  data.frame(run = "PSD5-only (2v2)",  event_type = psd5_run$summ$event_type,
             direction = "Inclusion UP", n = psd5_run$summ$n_up),
  data.frame(run = "PSD5-only (2v2)",  event_type = psd5_run$summ$event_type,
             direction = "Inclusion DOWN", n = psd5_run$summ$n_down))
p_hc <- ggplot(hc, aes(x = event_type, y = n, fill = direction)) +
  geom_col(position = "dodge", width = 0.8) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Inclusion UP" = "#d73027", "Inclusion DOWN" = "#2166ac")) +
  facet_wrap(~ run, nrow = 1) +
  labs(title = sprintf("rMATS hit counts (FDR<%.2f & |dPSI|>=%.2f)", FDR_THR, DPSI_THR),
       x = "Event type", y = "# sig events", fill = "") +
  theme_bw(base_size = 11)

p_ov <- ggplot(rbind(
  data.frame(event_type = overlap$event_type, category = "shared (both runs sig)", n = overlap$n_shared_sig),
  data.frame(event_type = overlap$event_type, category = "full-only sig",          n = overlap$n_full_sig - overlap$n_shared_sig),
  data.frame(event_type = overlap$event_type, category = "PSD5-only only sig",     n = overlap$n_PSD5_sig - overlap$n_shared_sig)
), aes(x = event_type, y = n, fill = category)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("shared (both runs sig)" = "#1f7a4d",
                                "full-only sig" = "#9ecae1",
                                "PSD5-only only sig" = "#fdae6b")) +
  labs(title = "Overlap of significant events: PSD5+8 full vs PSD5-only",
       x = "Event type", y = "# sig events", fill = "") +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_DIR, "rmats_summary_combined_PSD58.pdf"),
       (p_hc / p_ov) + plot_layout(heights = c(1, 1)),
       width = 11, height = 9)
cat("\nWrote rMATS summary plots + tables to:", OUT_DIR, "\n")
