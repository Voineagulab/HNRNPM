#!/usr/bin/env Rscript
# =============================================================================
# neuroblastoma_Fuchs2023.R
# Overlap of the 36 HNRNPM-regulated ("independent-backsplicing") circRNAs with
# neuroblastoma (NB) circRNA sets from Fuchs et al. 2023 (Nat Commun,
# doi:10.1038/s41467-023-38747-4). TWO disease-associated categories:
#   (1) NB-specific circRNAs        (MOESM8) — has coordinates -> A + B levels
#   (2) MYCN-amplified DE circRNAs  (MOESM6) — host gene only   -> B level only
# Matching: (A) exact hg38 back-splice junction (+/-1 nt, chr+strand agree);
#           (B) host gene symbol.
# Outputs -> RESULTS/neuroblastoma_Fuchs2023/
# Run on rna2 RStudio.
# =============================================================================
suppressPackageStartupMessages({ library(readxl); library(UpSetR); library(grid) })

## ---- paths ------------------------------------------------------------------
BASE   <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript"
W      <- file.path(BASE, "w0_circ_Disease")
S5BX   <- file.path(BASE, "Supplementary_Table", "Supplementary_Tables_JW_2.xlsx")
IVCLS  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ/BSJ_FSJ_classification_additive.tsv"
INDIR  <- file.path(W, "neuroblastoma")
OUT    <- file.path(W, "RESULTS", "neuroblastoma_Fuchs2023")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
f_moesm <- function(n) list.files(INDIR, pattern = sprintf("MOESM%d_ESM.xlsx$", n),
                                  recursive = TRUE, full.names = TRUE)[1]

## ---- helpers ----------------------------------------------------------------
nchr <- function(x) tolower(sub("^chr", "", as.character(x)))
up   <- function(x) toupper(trimws(as.character(x)))
key  <- function(c, st, s, e) sprintf("%s|%s|%d|%d", c, st, s, e)
mk_env <- function(c, st, s, e) { E <- new.env(hash = TRUE, size = length(c)*2L+8L)
  for (i in seq_along(c)) assign(key(c[i], st[i], s[i], e[i]), TRUE, envir = E); E }
# does query match any target in env, allowing +/-1 on each end
hit_env <- function(qc, qst, qs, qe, E) vapply(seq_along(qc), function(i){
  for (ds in c(0L,1L,-1L)) for (de in c(0L,1L,-1L))
    if (exists(key(qc[i], qst[i], qs[i]+ds, qe[i]+de), envir = E, inherits = FALSE)) return(TRUE)
  FALSE }, logical(1))
# Fisher row
frow <- function(category, level, universe, k, n, K, N, note = "") {
  if (is.na(N) || N <= 0) return(data.frame(category, match_level = level, universe,
    N = NA, K_category = NA, n_HnrnpM = NA, k_overlap = NA, expected = NA, odds_ratio = NA,
    OR_CI95_low = NA, OR_CI95_high = NA, p_one_sided = NA, p_two_sided = NA, note,
    stringsAsFactors = FALSE))
  tab <- matrix(c(k, K-k, n-k, N-K-(n-k)), 2, byrow = TRUE)
  f1 <- fisher.test(tab, alternative = "greater"); f2 <- fisher.test(tab)
  data.frame(category, match_level = level, universe, N = N, K_category = K, n_HnrnpM = n,
    k_overlap = k, expected = round(n*K/N, 3), odds_ratio = round(unname(f2$estimate), 2),
    OR_CI95_low = round(f2$conf.int[1], 2), OR_CI95_high = round(f2$conf.int[2], 2),
    p_one_sided = signif(f1$p.value, 4), p_two_sided = signif(f2$p.value, 4), note,
    stringsAsFactors = FALSE)
}
parse_circ <- function(v){ m <- regmatches(v, regexec("^([^:]+):([0-9]+)-([0-9]+):([+-])$", v))
  ok <- vapply(m, length, 1L) == 5
  data.frame(c = ifelse(ok, nchr(vapply(m, function(z) if(length(z)==5) z[2] else NA, "")), NA),
             s = suppressWarnings(as.integer(vapply(m, function(z) if(length(z)==5) z[3] else NA, ""))),
             e = suppressWarnings(as.integer(vapply(m, function(z) if(length(z)==5) z[4] else NA, ""))),
             st = ifelse(ok, vapply(m, function(z) if(length(z)==5) z[5] else NA, ""), NA),
             stringsAsFactors = FALSE) }

## ---- the 36 (S5B sheet, exact copy kept) + HnrnpM universe ------------------
s5b <- as.data.frame(read_excel(S5BX, sheet = "S5B_independent_circRNAs", skip = 1))
s5b <- s5b[!is.na(s5b$circRNA) & s5b$circRNA != "", ]
pc  <- parse_circ(s5b$circRNA); s5b$c <- pc$c; s5b$s <- pc$s; s5b$e <- pc$e; s5b$st <- pc$st
s5b$gene <- up(s5b$gene_symbol)

cls <- read.delim(IVCLS, stringsAsFactors = FALSE, check.names = FALSE)
pcu <- parse_circ(cls$circRNA); okc <- !is.na(pcu$s)
cls <- cls[okc, ]; cls$c <- pcu$c[okc]; cls$s <- pcu$s[okc]; cls$e <- pcu$e[okc]; cls$st <- pcu$st[okc]
cls$gene <- up(cls$gene_symbol)
hn_genes <- sort(unique(cls$gene[cls$gene != "" & !is.na(cls$gene)]))
g36      <- sort(unique(s5b$gene[s5b$gene != "" & !is.na(s5b$gene)]))

## ---- Fuchs tables -----------------------------------------------------------
# MOESM4: NB tumour circRNA landscape (all detected) -> universe + sensitivity
d4 <- as.data.frame(read_excel(f_moesm(4)))
m4 <- regmatches(d4$locus, regexec("(chr[0-9A-Za-z]+)\\(([+-])\\):([0-9]+)-([0-9]+)", d4$locus))
o4 <- vapply(m4, length, 1L) == 5
nb <- unique(data.frame(c = nchr(vapply(m4[o4], `[`, "", 2)), st = vapply(m4[o4], `[`, "", 3),
       s = as.integer(vapply(m4[o4], `[`, "", 4)), e = as.integer(vapply(m4[o4], `[`, "", 5)),
       gene = up(d4$gene_name[o4]), stringsAsFactors = FALSE))
nb_genes <- sort(unique(nb$gene[nb$gene != "" & !is.na(nb$gene)]))

# MOESM8: NB-specific circRNAs (coordinates) -> category (1)
d8 <- as.data.frame(read_excel(f_moesm(8))); id8 <- d8[[1]]
m8 <- regmatches(id8, regexec("_chr([0-9A-Za-z]+)([+-])([0-9]+)-([0-9]+)", id8)); o8 <- vapply(m8, length, 1L) == 5
nbspec <- data.frame(c = nchr(vapply(m8[o8], `[`, "", 2)), st = vapply(m8[o8], `[`, "", 3),
       s = as.integer(vapply(m8[o8], `[`, "", 4)), e = as.integer(vapply(m8[o8], `[`, "", 5)),
       gene = up(d8$gene_name[o8]), circbase = as.character(d8[["circ_id (circbase.org)"]][o8]),
       stringsAsFactors = FALSE)
nbspec_genes <- sort(unique(nbspec$gene[nbspec$gene != ""]))

# MOESM6: MYCN-amplified DE circRNAs, per host gene (up/down + log2FC) -> category (2)
myc <- do.call(rbind, lapply(c("MNA circs up","MNA circs down"), function(sh){
  d <- as.data.frame(read_excel(f_moesm(6), sheet = sh))
  data.frame(gene = up(d$gene_name), log2FC = suppressWarnings(as.numeric(d$log2FoldChange)),
             padj = suppressWarnings(as.numeric(d$padj)), stringsAsFactors = FALSE) }))
myc <- myc[!is.na(myc$gene) & myc$gene != "", ]
myc_genes <- sort(unique(myc$gene))
myc_l2 <- setNames(myc$log2FC[!duplicated(myc$gene)], myc$gene[!duplicated(myc$gene)])

## ---- coordinate matching flags for the 36 ----------------------------------
E_nb     <- mk_env(nb$c, nb$st, nb$s, nb$e)          # NB-detected landscape
E_nbspec <- mk_env(nbspec$c, nbspec$st, nbspec$s, nbspec$e)
s5b$in_nb     <- hit_env(s5b$c, s5b$st, s5b$s, s5b$e, E_nb)
s5b$in_nbspec <- hit_env(s5b$c, s5b$st, s5b$s, s5b$e, E_nbspec)
cls$in_nb     <- hit_env(cls$c, cls$st, cls$s, cls$e, E_nb)
cls$in_nbspec <- hit_env(cls$c, cls$st, cls$s, cls$e, E_nbspec)
nb$in_nbspec  <- hit_env(nb$c,  nb$st,  nb$s,  nb$e,  E_nbspec)

## ==== S5B CSV with appended A/B + direction columns (two categories) =========
NBspec_A <- s5b$in_nbspec
NBspec_B <- s5b$gene %in% nbspec_genes
MYCN_A   <- rep(FALSE, nrow(s5b))                     # MOESM6 has no coordinates
MYCN_B   <- s5b$gene %in% myc_genes
out <- s5b[, 1:12]
out[["NBspecific_A_coordinate_match"]] <- NBspec_A
out[["NBspecific_A_direction"]]        <- ifelse(NBspec_A, "up (NB-specific)", "")
out[["NBspecific_B_hostgene_match"]]   <- NBspec_B
out[["NBspecific_B_direction"]]        <- ifelse(NBspec_B, "up (NB-specific)", "")
out[["MYCN_A_coordinate_match"]]       <- MYCN_A
out[["MYCN_A_direction"]]              <- rep("NA (MOESM6 host-gene only, no coordinates)", nrow(s5b))
out[["MYCN_B_hostgene_match"]]         <- MYCN_B
out[["MYCN_B_direction"]] <- ifelse(MYCN_B, sprintf("%s (log2FC=%.2f)",
                                    ifelse(myc_l2[s5b$gene] > 0, "up", "down"), myc_l2[s5b$gene]), "")
out[["disease"]] <- "Neuroblastoma"
write.csv(out, file.path(OUT, "S5B_36circRNAs_Fuchs2023_annotated.csv"), row.names = FALSE)

## ==== Fisher (mirror z260710_circDE36_vs_NB) =================================
U2 <- intersect(hn_genes, nb_genes)
summ <- rbind(
  frow("NB-specific circRNAs", "coordinate", "assay-matched: HnrnpM-tested AND NB-detected (MOESM4)",
       sum(s5b$in_nbspec), sum(s5b$in_nb), sum(cls$in_nb & cls$in_nbspec), sum(cls$in_nb)),
  frow("NB-specific circRNAs", "coordinate", "sensitivity: NB circRNA landscape (MOESM4)",
       sum(s5b$in_nbspec), sum(s5b$in_nb), sum(nb$in_nbspec), nrow(nb)),
  frow("MYCN-amplified DE circRNAs", "coordinate", "(n/a)", NA, NA, NA, NA,
       note = "NOT ASSESSABLE: MOESM6 reports per host gene only, no coordinates."),
  frow("NB-specific circRNAs", "host-gene", "assay-matched: HnrnpM-tested AND NB host genes (MOESM4)",
       length(Reduce(intersect, list(g36, nbspec_genes, U2))), length(intersect(g36, U2)),
       length(intersect(nbspec_genes, U2)), length(U2)),
  frow("NB-specific circRNAs", "host-gene", "sensitivity: NB circRNA host genes (MOESM4)",
       length(Reduce(intersect, list(g36, nbspec_genes, nb_genes))), length(intersect(g36, nb_genes)),
       length(intersect(nbspec_genes, nb_genes)), length(nb_genes)),
  frow("MYCN-amplified DE circRNAs", "host-gene", "assay-matched: HnrnpM-tested AND NB host genes (MOESM4)",
       length(Reduce(intersect, list(g36, myc_genes, U2))), length(intersect(g36, U2)),
       length(intersect(myc_genes, U2)), length(U2)),
  frow("MYCN-amplified DE circRNAs", "host-gene", "sensitivity: NB circRNA host genes (MOESM4)",
       length(Reduce(intersect, list(g36, myc_genes, nb_genes))), length(intersect(g36, nb_genes)),
       length(intersect(myc_genes, nb_genes)), length(nb_genes)))
write.csv(summ, file.path(OUT, "fisher_summary.csv"), row.names = FALSE)

## ==== overlap tables ========================================================
ov_coord <- data.frame(host_gene = s5b$gene[s5b$in_nbspec],
  circRNA_hg38 = s5b$circRNA[s5b$in_nbspec], matched_set = "NB-specific (MOESM8)",
  stringsAsFactors = FALSE)
write.csv(ov_coord, file.path(OUT, "overlap_A_coordinate.csv"), row.names = FALSE)
gi_spec <- intersect(g36, nbspec_genes); gi_myc <- intersect(g36, myc_genes)
ov_gene <- rbind(
  if (length(gi_spec)) data.frame(host_gene = gi_spec, category = "NB-specific (MOESM8)",
             direction = "up (NB-specific)", stringsAsFactors = FALSE),
  if (length(gi_myc)) data.frame(host_gene = gi_myc, category = "MYCN-amplified DE (MOESM6)",
             direction = sprintf("%s (log2FC=%.2f)", ifelse(myc_l2[gi_myc] > 0, "up", "down"), myc_l2[gi_myc]),
             stringsAsFactors = FALSE))
write.csv(ov_gene, file.path(OUT, "overlap_B_hostgene.csv"), row.names = FALSE)

## ==== UpSet plots (PDF) =====================================================
# coordinate level: canonical anchoring so +/-1 variants align across sets
E_can <- new.env(hash = TRUE, size = nrow(cls)*20L)
for (i in seq_len(nrow(cls))) { kk <- key(cls$c[i], cls$st[i], cls$s[i], cls$e[i])
  for (ds in c(0L,1L,-1L)) for (de in c(0L,1L,-1L)) assign(key(cls$c[i], cls$st[i], cls$s[i]+ds, cls$e[i]+de), kk, envir = E_can) }
anchor <- function(c, st, s, e){ v <- key(c, st, s, e)
  if (exists(v, envir = E_can, inherits = FALSE)) get(v, envir = E_can) else v }
set36  <- unique(mapply(anchor, s5b$c, s5b$st, s5b$s, s5b$e))
setSpec<- unique(mapply(anchor, nbspec$c, nbspec$st, nbspec$s, nbspec$e))
setDet <- unique(mapply(anchor, nb$c, nb$st, nb$s, nb$e))
Lc <- list(`HnrnpM-36` = set36, `NB-specific` = setSpec, `NB-detected` = setDet)
pdf(file.path(OUT, "upset_coordinate.pdf"), width = 8, height = 5)
print(upset(fromList(Lc), sets = rev(names(Lc)), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Shared circRNAs (exact junction)", sets.x.label = "circRNAs per set", text.scale = 1.2))
dev.off()
Lg <- list(`HnrnpM-36` = g36, `NB-specific` = nbspec_genes, `MYCN-DE` = myc_genes)
pdf(file.path(OUT, "upset_hostgene.pdf"), width = 8, height = 5)
print(upset(fromList(Lg), sets = rev(names(Lg)), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Shared host genes", sets.x.label = "genes per set", text.scale = 1.2))
dev.off()

cat(sprintf("Fuchs done. NB-specific coord hits=%d; NB-specific gene hits=%d; MYCN gene hits=%d\n",
            sum(NBspec_A), length(gi_spec), length(gi_myc)))
