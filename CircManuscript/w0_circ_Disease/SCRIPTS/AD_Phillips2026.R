#!/usr/bin/env Rscript
# =============================================================================
# AD_Phillips2026.R
# Overlap of the 36 HNRNPM-regulated circRNAs with blood AD-associated circRNAs
# from Phillips et al. 2026 (blood-based circRNA AD diagnosis).
# Coordinates ARE reported (Location = chr_start_end, hg38) -> matched at BOTH
# the coordinate (A) and host-gene (B) level.
# DISEASE-ASSOCIATED = circRNAs significant at FDR<0.05 in the AD-association
# analyses (Supplementary Table S3 = DCC, S4 = CIRI2). Direction from log2FC
# (up = enriched in AD blood). "Merely detected but not significant" circRNAs
# are used ONLY to build the Fisher universe, never flagged A/B.
# Outputs -> RESULTS/AD_Phillips2026/
# Run on rna2 RStudio.
# =============================================================================
suppressPackageStartupMessages({ library(readxl); library(UpSetR) })

BASE  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript"
W     <- file.path(BASE, "w0_circ_Disease")
S5BX  <- file.path(BASE, "Supplementary_Table", "Supplementary_Tables_JW_2.xlsx")
IVCLS <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ/BSJ_FSJ_classification_additive.tsv"
XLSX  <- list.files(file.path(W, "Alzheimer", "Phillips2026"), pattern = "\\.xlsx$", full.names = TRUE)
XLSX  <- XLSX[!grepl("/\\._", XLSX)][1]
OUT   <- file.path(W, "RESULTS", "AD_Phillips2026")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

nchr <- function(x) tolower(sub("^chr", "", as.character(x)))
up   <- function(x) toupper(trimws(as.character(x)))
frow <- function(category, level, universe, k, n, K, N, note = "") {
  if (is.na(N) || N <= 0 || (N-K-(n-k)) < 0) return(data.frame(category, match_level = level, universe, N = N,
    K_category = K, n_HnrnpM = n, k_overlap = k, expected = NA, odds_ratio = NA, OR_CI95_low = NA,
    OR_CI95_high = NA, p_one_sided = NA, p_two_sided = NA, note = paste(note,"invalid"), stringsAsFactors = FALSE))
  tab <- matrix(c(k, K-k, n-k, N-K-(n-k)), 2, byrow = TRUE)
  f1 <- fisher.test(tab, alternative = "greater"); f2 <- fisher.test(tab)
  data.frame(category, match_level = level, universe, N = N, K_category = K, n_HnrnpM = n, k_overlap = k,
    expected = round(n*K/N, 3), odds_ratio = round(unname(f2$estimate), 2),
    OR_CI95_low = round(f2$conf.int[1], 2), OR_CI95_high = round(f2$conf.int[2], 2),
    p_one_sided = signif(f1$p.value, 4), p_two_sided = signif(f2$p.value, 4), note, stringsAsFactors = FALSE)
}
parse_circ <- function(v){ m <- regmatches(v, regexec("^([^:]+):([0-9]+)-([0-9]+):([+-])$", v))
  ok <- vapply(m, length, 1L) == 5
  list(c = ifelse(ok, nchr(vapply(m, function(z) if(length(z)==5) z[2] else NA, "")), NA),
       s = suppressWarnings(as.integer(vapply(m, function(z) if(length(z)==5) z[3] else NA, ""))),
       e = suppressWarnings(as.integer(vapply(m, function(z) if(length(z)==5) z[4] else NA, ""))),
       st = ifelse(ok, vapply(m, function(z) if(length(z)==5) z[5] else NA, ""), NA)) }

## ---- the 36 + HnrnpM universe ----------------------------------------------
s5b <- as.data.frame(read_excel(S5BX, sheet = "S5B_independent_circRNAs", skip = 1))
s5b <- s5b[!is.na(s5b$circRNA) & s5b$circRNA != "", ]
p <- parse_circ(s5b$circRNA); s5b$c <- p$c; s5b$s <- p$s; s5b$e <- p$e; s5b$gene <- up(s5b$gene_symbol)
cls <- read.delim(IVCLS, stringsAsFactors = FALSE, check.names = FALSE)
pu <- parse_circ(cls$circRNA); ok <- !is.na(pu$s)
cls <- cls[ok, ]; cls$c <- pu$c[ok]; cls$s <- pu$s[ok]; cls$e <- pu$e[ok]; cls$gene <- up(cls$gene_symbol)
hn_genes <- sort(unique(cls$gene[cls$gene != ""])); g36 <- sort(unique(s5b$gene[s5b$gene != ""]))

## ---- Phillips association tables (S3 = DCC, S4 = CIRI2) ---------------------
read_ph <- function(sheet){
  d <- as.data.frame(read_excel(XLSX, sheet = sheet, skip = 1)); d <- d[!is.na(d$Location), ]
  loc <- strsplit(as.character(d$Location), "_")
  data.frame(gene = up(d$Gene),
             c = nchr(vapply(loc, `[`, "", 1)),
             s = suppressWarnings(as.integer(vapply(loc, function(z) z[2], ""))),
             e = suppressWarnings(as.integer(vapply(loc, function(z) z[3], ""))),
             log2FC = suppressWarnings(as.numeric(d$log2FC)),
             FDR = suppressWarnings(as.numeric(d$FDR)), stringsAsFactors = FALSE) }
s3 <- read_ph("Table S3"); s4 <- read_ph("Table S4")
detected <- unique(rbind(s3, s4)[, c("c","s","e","gene")])                 # tested/detected -> universe
assoc    <- unique(rbind(s3[s3$FDR < 0.05 & !is.na(s3$FDR), ], s4[s4$FDR < 0.05 & !is.na(s4$FDR), ]))
assoc    <- assoc[order(assoc$FDR), ]
assoc_genes <- sort(unique(assoc$gene[assoc$gene != ""]))
det_genes   <- sort(unique(detected$gene[detected$gene != ""]))

## ---- coordinate matching (chr + start/end +/-1; Phillips has no strand) -----
mk <- function(c, s, e){ E <- new.env(hash = TRUE, size = length(c)*2L+8L)
  for (i in seq_along(c)) assign(sprintf("%s|%d|%d", c[i], s[i], e[i]), TRUE, E); E }
hit <- function(qc, qs, qe, E) vapply(seq_along(qc), function(i){
  for (ds in c(0L,1L,-1L)) for (de in c(0L,1L,-1L))
    if (exists(sprintf("%s|%d|%d", qc[i], qs[i]+ds, qe[i]+de), E, inherits = FALSE)) return(TRUE)
  FALSE }, logical(1))
E_det <- mk(detected$c, detected$s, detected$e); E_as <- mk(assoc$c, assoc$s, assoc$e)
s5b$in_det <- hit(s5b$c, s5b$s, s5b$e, E_det); s5b$in_as <- hit(s5b$c, s5b$s, s5b$e, E_as)
cls$in_det <- hit(cls$c, cls$s, cls$e, E_det); cls$in_as <- hit(cls$c, cls$s, cls$e, E_as)
det_is_as  <- hit(detected$c, detected$s, detected$e, E_as)

# per-match direction (A: by coordinate; B: by gene, most significant)
Astat <- function(qc, qs, qe) vapply(seq_along(qc), function(i){
  for (ds in c(0L,1L,-1L)) for (de in c(0L,1L,-1L)){
    h <- which(assoc$c == qc[i] & assoc$s == qs[i]+ds & assoc$e == qe[i]+de)
    if (length(h)) return(sprintf("%s in AD blood (log2FC=%.2f, FDR=%.2g)",
                                  ifelse(assoc$log2FC[h[1]] > 0, "up", "down"), assoc$log2FC[h[1]], assoc$FDR[h[1]])) }
  "" }, "")
Bstat <- function(gene) vapply(gene, function(g){ h <- which(assoc$gene == g); if (!length(h)) return("")
  h <- h[order(assoc$FDR[h])][1]
  sprintf("%s in AD blood (log2FC=%.2f, FDR=%.2g)", ifelse(assoc$log2FC[h] > 0, "up", "down"),
          assoc$log2FC[h], assoc$FDR[h]) }, "")

## ==== S5B CSV with appended A/B + direction =================================
A <- s5b$in_as; B <- s5b$gene %in% assoc_genes
out <- s5b[, 1:12]
out[["AD_A_coordinate_match"]] <- A
out[["AD_A_direction"]]        <- ifelse(A, Astat(s5b$c, s5b$s, s5b$e), "")
out[["AD_B_hostgene_match"]]   <- B
out[["AD_B_direction"]]        <- ifelse(B, Bstat(s5b$gene), "")
out[["disease"]] <- "Alzheimer"
write.csv(out, file.path(OUT, "S5B_36circRNAs_Phillips2026_annotated.csv"), row.names = FALSE)

## ==== Fisher ================================================================
U2 <- intersect(hn_genes, det_genes)
summ <- rbind(
  frow("AD-associated blood circRNAs (FDR<0.05)", "coordinate",
       "assay-matched: HnrnpM-tested AND Phillips-detected",
       sum(s5b$in_as), sum(s5b$in_det), sum(cls$in_det & cls$in_as), sum(cls$in_det)),
  frow("AD-associated blood circRNAs (FDR<0.05)", "coordinate",
       "sensitivity: Phillips detected circRNAs (S3 union S4)",
       sum(s5b$in_as), sum(s5b$in_det), sum(det_is_as), nrow(detected)),
  frow("AD-associated blood circRNAs (FDR<0.05)", "host-gene",
       "assay-matched: HnrnpM-tested AND Phillips-detected host genes",
       length(Reduce(intersect, list(g36, assoc_genes, U2))), length(intersect(g36, U2)),
       length(intersect(assoc_genes, U2)), length(U2)),
  frow("AD-associated blood circRNAs (FDR<0.05)", "host-gene",
       "sensitivity: Phillips-detected host genes",
       length(Reduce(intersect, list(g36, assoc_genes, det_genes))), length(intersect(g36, det_genes)),
       length(intersect(assoc_genes, det_genes)), length(det_genes)))
write.csv(summ, file.path(OUT, "fisher_summary.csv"), row.names = FALSE)

## ==== overlap tables ========================================================
oc <- s5b[s5b$in_as, ]
write.csv(data.frame(host_gene = oc$gene, circRNA_hg38 = oc$circRNA,
                     AD_direction = Astat(oc$c, oc$s, oc$e), stringsAsFactors = FALSE),
          file.path(OUT, "overlap_A_coordinate.csv"), row.names = FALSE)
gi <- intersect(g36, assoc_genes)
write.csv(data.frame(host_gene = gi, AD_direction = Bstat(gi), stringsAsFactors = FALSE),
          file.path(OUT, "overlap_B_hostgene.csv"), row.names = FALSE)

## ==== UpSet (coordinate, PDF) ===============================================
kf <- function(c, s, e) sprintf("%s|%d|%d", c, s, e)
E_can <- new.env(hash = TRUE, size = nrow(cls)*20L)
for (i in seq_len(nrow(cls))) { kk <- kf(cls$c[i], cls$s[i], cls$e[i])
  for (ds in c(0L,1L,-1L)) for (de in c(0L,1L,-1L)) assign(kf(cls$c[i], cls$s[i]+ds, cls$e[i]+de), kk, E_can) }
anc <- function(c, s, e){ v <- kf(c, s, e); if (exists(v, E_can, inherits = FALSE)) get(v, E_can) else v }
L <- list(`HnrnpM-36` = unique(mapply(anc, s5b$c, s5b$s, s5b$e)),
          `AD-associated` = unique(mapply(anc, assoc$c, assoc$s, assoc$e)),
          `AD-detected` = unique(mapply(anc, detected$c, detected$s, detected$e)))
pdf(file.path(OUT, "upset_coordinate.pdf"), width = 8, height = 5)
print(upset(fromList(L), sets = rev(names(L)), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Shared circRNAs (exact junction)", sets.x.label = "circRNAs per set", text.scale = 1.2))
dev.off()

cat(sprintf("Phillips done. detected=%d assoc(FDR<0.05)=%d ; coord A hits=%d ; host-gene B hits=%d (%s)\n",
            nrow(detected), nrow(assoc), sum(A), length(gi), paste(gi, collapse = ", ")))
