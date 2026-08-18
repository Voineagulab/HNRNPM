#!/usr/bin/env Rscript
# =============================================================================
# glioblastoma_Song2016.R
# Overlap of the 36 HNRNPM-regulated circRNAs with glioma circRNAs from
# Song et al. 2016 (NAR, UROBORUS).
#
# DISEASE-ASSOCIATION IS RE-DERIVED FROM THE PRIMARY DATA (see README):
#   The paper does NOT publish per-circRNA fold changes. Supplementary Table S5
#   gives all 572 highly-expressed circRNAs (hg19: chr/start/end/gene) plus the
#   per-sample RPM matrix (RL1..RL46); Supplementary Table S4 assigns samples to
#   groups (20 glioblastoma [GBM], 19 normal, 7 oligodendroglioma). We reproduce
#   the authors' own test: Wilcoxon rank-sum (asymptotic; == scipy.stats.ranksums)
#   GBM vs normal on the RPM values, Benjamini-Hochberg correction. A circRNA is
#   "disease-associated" if q < 0.05. Direction (enriched/depleted in GBM) = sign
#   of (mean GBM RPM - mean normal RPM). Since no FC is published, we report a
#   DERIVED effect size log2((meanGBM+eps)/(meanNorm+eps)) as a magnitude
#   indicator only (clearly labelled "derived, not published").
# Coordinates are lifted hg19 -> hg38 (rtracklayer) before matching the 36.
# Outputs -> RESULTS/glioblastoma_Song2016/
# Run on rna2 RStudio (needs internet for the UCSC chain).
# =============================================================================
suppressPackageStartupMessages({ library(readxl); library(UpSetR); library(rtracklayer); library(GenomicRanges) })

BASE  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript"
W     <- file.path(BASE, "w0_circ_Disease")
S5BX  <- file.path(BASE, "Supplementary_Table", "Supplementary_Tables_JW_2.xlsx")
IVCLS <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ/BSJ_FSJ_classification_additive.tsv"
XLSX  <- list.files(file.path(W, "glioblastoma"), pattern = "UROBORUS_SupTables.xlsx$",
                    recursive = TRUE, full.names = TRUE)[1]
OUT   <- file.path(W, "RESULTS", "glioblastoma_Song2016")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

nchr <- function(x) tolower(sub("^chr", "", as.character(x)))
up   <- function(x) toupper(trimws(as.character(x)))
key  <- function(c, st, s, e) sprintf("%s|%s|%d|%d", c, st, s, e)
parse_circ <- function(v){ m <- regmatches(v, regexec("^([^:]+):([0-9]+)-([0-9]+):([+-])$", v))
  ok <- vapply(m, length, 1L) == 5
  list(c = ifelse(ok, nchr(vapply(m, function(z) if(length(z)==5) z[2] else NA, "")), NA),
       s = suppressWarnings(as.integer(vapply(m, function(z) if(length(z)==5) z[3] else NA, ""))),
       e = suppressWarnings(as.integer(vapply(m, function(z) if(length(z)==5) z[4] else NA, ""))),
       st = ifelse(ok, vapply(m, function(z) if(length(z)==5) z[5] else NA, ""), NA)) }
frow <- function(category, level, universe, k, n, K, N, note = "") {
  if (is.na(N) || N <= 0) return(data.frame(category, match_level = level, universe, N = NA,
    K_category = NA, n_HnrnpM = NA, k_overlap = NA, expected = NA, odds_ratio = NA, OR_CI95_low = NA,
    OR_CI95_high = NA, p_one_sided = NA, p_two_sided = NA, note, stringsAsFactors = FALSE))
  tab <- matrix(c(k, K-k, n-k, N-K-(n-k)), 2, byrow = TRUE)
  f1 <- fisher.test(tab, alternative = "greater"); f2 <- fisher.test(tab)
  data.frame(category, match_level = level, universe, N = N, K_category = K, n_HnrnpM = n, k_overlap = k,
    expected = round(n*K/N, 3), odds_ratio = round(unname(f2$estimate), 2),
    OR_CI95_low = round(f2$conf.int[1], 2), OR_CI95_high = round(f2$conf.int[2], 2),
    p_one_sided = signif(f1$p.value, 4), p_two_sided = signif(f2$p.value, 4), note, stringsAsFactors = FALSE)
}
# scipy.stats.ranksums equivalent (asymptotic; average ranks for ties; no tie/continuity corr.)
ranksums_p <- function(x, y){ x <- x[!is.na(x)]; y <- y[!is.na(y)]
  n1 <- length(x); n2 <- length(y); N <- n1 + n2; r <- rank(c(x, y)); s <- sum(r[seq_len(n1)])
  2 * pnorm(-abs((s - n1*(N+1)/2) / sqrt(n1*n2*(N+1)/12))) }

## ---- the 36 + HnrnpM universe (hg38) ---------------------------------------
s5b <- as.data.frame(read_excel(S5BX, sheet = "S5B_independent_circRNAs", skip = 1))
s5b <- s5b[!is.na(s5b$circRNA) & s5b$circRNA != "", ]
p <- parse_circ(s5b$circRNA); s5b$c <- p$c; s5b$s <- p$s; s5b$e <- p$e; s5b$st <- p$st; s5b$gene <- up(s5b$gene_symbol)
cls <- read.delim(IVCLS, stringsAsFactors = FALSE, check.names = FALSE)
pu <- parse_circ(cls$circRNA); ok <- !is.na(pu$s)
cls <- cls[ok, ]; cls$c <- pu$c[ok]; cls$s <- pu$s[ok]; cls$e <- pu$e[ok]; cls$st <- pu$st[ok]; cls$gene <- up(cls$gene_symbol)
hn_genes <- sort(unique(cls$gene[cls$gene != ""])); g36 <- sort(unique(s5b$gene[s5b$gene != ""]))

## ---- Song DE: Wilcoxon rank-sum + BH, GBM vs normal (from S5 + S4) ----------
S5 <- as.data.frame(read_excel(XLSX, sheet = "Table S5"))
S4 <- as.data.frame(read_excel(XLSX, sheet = "Table S4")); names(S4)[1:2] <- c("sid","tissue")
S4$sid <- trimws(S4$sid); S4$tissue <- tolower(trimws(S4$tissue))
samp <- names(S5)[5:ncol(S5)]; grp <- setNames(S4$tissue, S4$sid)[samp]
gbm <- samp[grepl("^glio", grp)]; nrm <- samp[grp == "normal"]
X <- sapply(S5[, gbm], as.numeric); Y <- sapply(S5[, nrm], as.numeric)
pv <- vapply(seq_len(nrow(S5)), function(i) ranksums_p(X[i, ], Y[i, ]), numeric(1))
md <- rowMeans(X, na.rm = TRUE) - rowMeans(Y, na.rm = TRUE)
q  <- p.adjust(pv, method = "BH")
eps <- 1e-3
l2fc <- log2((rowMeans(X, na.rm = TRUE) + eps) / (rowMeans(Y, na.rm = TRUE) + eps))   # DERIVED, not published
de <- data.frame(chr = nchr(S5$chr), s = as.integer(S5$start), e = as.integer(S5$end),
                 gene = up(S5$gene), q = q, l2fc = l2fc, dir = ifelse(md > 0, "up in GBM", "down in GBM"),
                 stringsAsFactors = FALSE)
cat(sprintf("Song DE (GBM %d vs normal %d): q<0.05 down=%d up=%d ; of 572\n",
            length(gbm), length(nrm), sum(q < 0.05 & md < 0), sum(q < 0.05 & md > 0)))

## ---- liftOver hg19 -> hg38 (all 572; carry DE stats) -----------------------
chain_gz <- file.path(OUT, "hg19ToHg38.over.chain.gz"); chain_f <- file.path(OUT, "hg19ToHg38.over.chain")
if (!file.exists(chain_f)) {
  if (!file.exists(chain_gz)) download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", chain_gz, quiet = TRUE)
  R.utils::gunzip(chain_gz, destname = chain_f, remove = FALSE, overwrite = TRUE) }
ch <- import.chain(chain_f)
gr <- GRanges(paste0("chr", de$chr), IRanges(de$s + 1, de$e)); mcols(gr)$idx <- seq_len(nrow(de))   # 0-based -> 1-based
lo <- liftOver(gr, ch); mapped <- lengths(lo) > 0; lo1 <- unlist(range(lo))
cat(sprintf("liftOver: %d/%d mapped to hg38\n", sum(mapped), length(gr)))
gl <- de[mapped, ]; gl$c <- nchr(as.character(seqnames(lo1))); gl$s38 <- start(lo1) - 1; gl$e38 <- end(lo1)   # back to 0-based
gl$st <- s5b$st[match(gl$gene, s5b$gene)]   # strand unknown in S5; only used for keying (NA tolerated)
glDE <- gl[gl$q < 0.05, ]                    # disease-associated set (hg38)

## ---- coordinate matching (chr + start/end +/-1; strand ignored, S5 has none)
mk_env2 <- function(c, s, e){ E <- new.env(hash = TRUE, size = length(c)*2L+8L)
  for (i in seq_along(c)) assign(sprintf("%s|%d|%d", c[i], s[i], e[i]), TRUE, E); E }
hit2 <- function(qc, qs, qe, E) vapply(seq_along(qc), function(i){
  for (ds in c(0L,1L,-1L)) for (de_ in c(0L,1L,-1L))
    if (exists(sprintf("%s|%d|%d", qc[i], qs[i]+ds, qe[i]+de_), E, inherits = FALSE)) return(TRUE)
  FALSE }, logical(1))
E_det <- mk_env2(gl$c, gl$s38, gl$e38)       # all glioma-detected (hg38)
E_de  <- mk_env2(glDE$c, glDE$s38, glDE$e38) # glioma disease-associated (q<0.05)
s5b$in_det <- hit2(s5b$c, s5b$s, s5b$e, E_det); s5b$in_de <- hit2(s5b$c, s5b$s, s5b$e, E_de)
cls$in_det <- hit2(cls$c, cls$s, cls$e, E_det); cls$in_de <- hit2(cls$c, cls$s, cls$e, E_de)
gl$is_de   <- gl$q < 0.05
de_genes  <- sort(unique(glDE$gene[glDE$gene != ""]))
det_genes <- sort(unique(gl$gene[gl$gene != ""]))

# per-matched direction/FC for the 36 (nearest DE circRNA by coordinate; and by gene for B)
match_stat <- function(qc, qs, qe){ vapply(seq_along(qc), function(i){
  for (ds in c(0L,1L,-1L)) for (de_ in c(0L,1L,-1L)){
    h <- which(glDE$c == qc[i] & glDE$s38 == qs[i]+ds & glDE$e38 == qe[i]+de_)
    if (length(h)) return(sprintf("%s (log2FC=%.2f, q=%.2g; derived)", glDE$dir[h[1]], glDE$l2fc[h[1]], glDE$q[h[1]])) }
  "" }, "") }
Bstat <- function(gene){ vapply(gene, function(g){
  h <- which(glDE$gene == g); if (!length(h)) return("")
  h <- h[order(glDE$q[h])][1]
  sprintf("%s (log2FC=%.2f, q=%.2g; derived)", glDE$dir[h], glDE$l2fc[h], glDE$q[h]) }, "") }

## ==== S5B CSV with appended A/B + direction ================================
A <- s5b$in_de; B <- s5b$gene %in% de_genes
out <- s5b[, 1:12]
out[["glioma_A_coordinate_match"]] <- A
out[["glioma_A_direction"]]        <- ifelse(A, match_stat(s5b$c, s5b$s, s5b$e), "")
out[["glioma_B_hostgene_match"]]   <- B
out[["glioma_B_direction"]]        <- ifelse(B, Bstat(s5b$gene), "")
out[["disease"]] <- "Glioblastoma"
write.csv(out, file.path(OUT, "S5B_36circRNAs_Song2016_annotated.csv"), row.names = FALSE)

## ==== Fisher ================================================================
U2 <- intersect(hn_genes, det_genes)
summ <- rbind(
  frow("Glioma DE circRNAs (q<0.05, GBM vs normal)", "coordinate",
       "assay-matched: HnrnpM-tested AND glioma-detected (S5, lifted)",
       sum(s5b$in_de), sum(s5b$in_det), sum(cls$in_det & cls$in_de), sum(cls$in_det)),
  frow("Glioma DE circRNAs (q<0.05, GBM vs normal)", "coordinate",
       "sensitivity: glioma highly-expressed landscape (572, S5)",
       sum(s5b$in_de), sum(s5b$in_det), sum(gl$is_de), nrow(gl)),
  frow("Glioma DE circRNAs (q<0.05, GBM vs normal)", "host-gene",
       "assay-matched: HnrnpM-tested AND glioma-detected host genes",
       length(Reduce(intersect, list(g36, de_genes, U2))), length(intersect(g36, U2)),
       length(intersect(de_genes, U2)), length(U2)),
  frow("Glioma DE circRNAs (q<0.05, GBM vs normal)", "host-gene",
       "sensitivity: glioma-detected host genes (S5)",
       length(Reduce(intersect, list(g36, de_genes, det_genes))), length(intersect(g36, det_genes)),
       length(intersect(de_genes, det_genes)), length(det_genes)))
write.csv(summ, file.path(OUT, "fisher_summary.csv"), row.names = FALSE)

## ==== overlap tables ========================================================
oc <- s5b[s5b$in_de, ]
write.csv(data.frame(host_gene = oc$gene, circRNA_hg38 = oc$circRNA,
                     glioma_direction = match_stat(oc$c, oc$s, oc$e), stringsAsFactors = FALSE),
          file.path(OUT, "overlap_A_coordinate.csv"), row.names = FALSE)
gi <- intersect(g36, de_genes)
write.csv(data.frame(host_gene = gi, glioma_direction = Bstat(gi), stringsAsFactors = FALSE),
          file.path(OUT, "overlap_B_hostgene.csv"), row.names = FALSE)

## ==== UpSet (PDF) ===========================================================
kf <- function(c, s, e) sprintf("%s|%d|%d", c, s, e)
E_can <- new.env(hash = TRUE, size = nrow(cls)*20L)
for (i in seq_len(nrow(cls))) { kk <- kf(cls$c[i], cls$s[i], cls$e[i])
  for (ds in c(0L,1L,-1L)) for (de_ in c(0L,1L,-1L)) assign(kf(cls$c[i], cls$s[i]+ds, cls$e[i]+de_), kk, E_can) }
anc <- function(c, s, e){ v <- kf(c, s, e); if (exists(v, E_can, inherits = FALSE)) get(v, E_can) else v }
L <- list(`HnrnpM-36` = unique(mapply(anc, s5b$c, s5b$s, s5b$e)),
          `glioma-DE`  = unique(mapply(anc, glDE$c, glDE$s38, glDE$e38)),
          `glioma-detected` = unique(mapply(anc, gl$c, gl$s38, gl$e38)))
pdf(file.path(OUT, "upset_coordinate.pdf"), width = 8, height = 5)
print(upset(fromList(L), sets = rev(names(L)), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Shared circRNAs (exact junction)", sets.x.label = "circRNAs per set", text.scale = 1.2))
dev.off()

cat(sprintf("Song done. coord A hits=%d ; host-gene B hits=%d\n", sum(A), length(gi)))
