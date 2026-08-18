#!/usr/bin/env Rscript
# =============================================================================
# neuroblastoma_Fuchs2023_MYNC.R
# HOST-GENE-LEVEL test of whether the 36 LOI-circRNA host genes are over-
# represented among the host genes whose circRNA abundance is REDUCED in
# MYCN-amplified (MNA) neuroblastoma, per Fuchs et al. 2023 (Nat Commun,
# doi:10.1038/s41467-023-38747-4), MOESM6 sheet "MNA circs down".
#
# DESIGN ("down, both split")
#   * category   : MOESM6 "MNA circs down" only. The "up" sheet is read for
#                  annotation but is NOT part of the Fisher test.
#   * identifier : Ensembl gene ID throughout. No symbol -> ENSG conversion is
#                  performed or needed: our own pipeline output carries a
#                  gene_id column. Version suffixes are stripped.
#   * SPLIT      : some circRNAs are annotated to more than one gene, and carry
#                  a comma-separated gene_id (e.g. BBS9 =
#                  "ENSG00000122507,ENSG00000238090"). BOTH the background and
#                  the 36 are split on those separators, so each Ensembl ID is
#                  counted once. Splitting one side only would bias the odds
#                  ratio upward.
#   * universe   : the host gene IDs of all 5,820 HnrnpM-tested circRNAs. There
#                  is no additional requirement that a gene also be detected in
#                  the Fuchs tumour cohort (MOESM4 is not used).
#   Fuchs needs no splitting: MOESM6 carries exactly one Ensembl ID per row.
#
# EXPECTED (the script stops if these are not reproduced)
#   N = 3025 ; K = 296 ; n = 35 ; k = 8 ; OR = 2.78 (95% CI 1.08-6.36) ; P = 0.017
#
# Outputs -> RESULTS/neuroblastoma_Fuchs2023/MYNC/
# Run on rna2 RStudio.
# =============================================================================
suppressPackageStartupMessages({ library(readxl); library(UpSetR); library(grid) })

## ---- paths ------------------------------------------------------------------
BASE   <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript"
W      <- file.path(BASE, "w0_circ_Disease")
S5BX   <- file.path(BASE, "Supplementary_Table", "Supplementary_Tables_JW_2.xlsx")
# background + gene_id source. NB: BSJ_FSJ_classification_additive.tsv carries the
# same 5,820 circRNAs but has NO gene_id column, which is why the BSJ limma table
# is read instead.
IVBSJ  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv"
INDIR  <- file.path(W, "neuroblastoma")
OUT    <- file.path(W, "RESULTS", "neuroblastoma_Fuchs2023", "MYNC")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
# skip Excel lock files (~$...) and AppleDouble (._...) siblings, which also end
# in "MOESM<n>_ESM.xlsx"
f_moesm <- function(n){
  f <- list.files(INDIR, pattern = sprintf("MOESM%d_ESM\\.xlsx$", n),
                  recursive = TRUE, full.names = TRUE)
  f <- f[!grepl("^[~.]", basename(f))]
  if (!length(f)) stop("MOESM", n, " workbook not found under ", INDIR)
  f[1] }

## ---- helpers ----------------------------------------------------------------
sv  <- function(x) sub("\\.[0-9]+$", "", trimws(as.character(x)))   # drop ENSG version
# split a gene_id cell into its individual Ensembl IDs (THE "split" step)
ids <- function(x){ p <- strsplit(trimws(as.character(x)), "[,;|/[:space:]]+")
  lapply(p, function(v){ v <- sv(v); unique(v[grepl("^ENSG[0-9]+$", v)]) }) }
# Fisher row (identical to neuroblastoma_Fuchs2023.R)
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

## ---- the 36 (S5B sheet, exact copy kept) ------------------------------------
s5b <- as.data.frame(read_excel(S5BX, sheet = "S5B_independent_circRNAs", skip = 1))
s5b <- s5b[!is.na(s5b$circRNA) & s5b$circRNA != "", ]

## ---- background: host gene IDs of all HnrnpM-tested circRNAs ---------------
bsj <- read.delim(IVBSJ, stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene_id" %in% names(bsj)) stop("no gene_id column in ", IVBSJ)
bsj_ids <- ids(bsj$gene_id)                       # list, one character vector per circRNA
names(bsj_ids) <- bsj$circRNA
BG <- sort(unique(unlist(bsj_ids, use.names = FALSE)))          # N

## ---- the 36's host gene IDs (same split) ------------------------------------
miss <- setdiff(s5b$circRNA, names(bsj_ids))
if (length(miss)) stop("circRNAs absent from the BSJ table: ", paste(miss, collapse = ", "))
s5b_ids   <- bsj_ids[s5b$circRNA]
s5b$gene_ids   <- vapply(s5b_ids, paste, "", collapse = ",")
s5b$n_gene_ids <- vapply(s5b_ids, length, 1L)
G36 <- sort(unique(unlist(s5b_ids, use.names = FALSE)))         # n (before universe cap)

## ---- Fuchs MOESM6: MNA circs down (category) + up (annotation only) --------
rd6 <- function(sh){ d <- as.data.frame(read_excel(f_moesm(6), sheet = sh))
  data.frame(ens = sv(d[["ensembl ID"]]), fuchs_gene_name = as.character(d$gene_name),
             log2FC = suppressWarnings(as.numeric(d$log2FoldChange)),
             padj = suppressWarnings(as.numeric(d$padj)), fuchs_set = sh,
             stringsAsFactors = FALSE) }
dn <- rd6("MNA circs down"); upx <- rd6("MNA circs up")
mo <- rbind(dn, upx)
mo <- mo[!is.na(mo$ens) & nzchar(mo$ens), ]
mo <- mo[order(mo$padj), ]; mo <- mo[!duplicated(mo$ens), ]   # most significant per gene
DOWN <- sort(unique(dn$ens[nzchar(dn$ens)]))                  # THE tested category
UPS  <- sort(unique(upx$ens[nzchar(upx$ens)]))

## ==== Fisher: down-regulated only, both sides split =========================
N <- length(BG)
K <- length(intersect(DOWN, BG))
n <- length(intersect(G36,  BG))
k <- length(Reduce(intersect, list(G36, DOWN, BG)))
summ <- frow("MYCN-amplified DOWN-regulated circRNA host genes (MOESM6 'MNA circs down')",
             "host-gene (Ensembl gene ID, multi-gene assignments split)",
             "HnrnpM-tested circRNA host genes (no MOESM4 detectability filter)",
             k, n, K, N)
write.csv(summ, file.path(OUT, "fisher_summary.csv"), row.names = FALSE)

## ---- expected values: fail loudly if the design drifts ---------------------
EXP <- list(N = 3025L, K = 296L, n = 35L, k = 8L, OR = 2.78, LO = 1.08, HI = 6.36, P = 0.0171)
chk <- c(N = N == EXP$N, K = K == EXP$K, n = n == EXP$n, k = k == EXP$k,
         OR = abs(summ$odds_ratio - EXP$OR) < 0.005,
         CI_low = abs(summ$OR_CI95_low - EXP$LO) < 0.005,
         CI_high = abs(summ$OR_CI95_high - EXP$HI) < 0.005,
         P = abs(summ$p_one_sided - EXP$P) < 0.0005)
if (!all(chk)) warning("DESIGN CHECK FAILED for: ", paste(names(chk)[!chk], collapse = ", "),
  "\n  got   N=", N, " K=", K, " n=", n, " k=", k, " OR=", summ$odds_ratio,
  " CI=", summ$OR_CI95_low, "-", summ$OR_CI95_high, " P=", summ$p_one_sided,
  "\n  expect N=3025 K=296 n=35 k=8 OR=2.78 CI=1.08-6.36 P=0.0171")

## ==== S5B CSV with appended match columns ===================================
hit <- function(v, S) vapply(s5b_ids, function(g) any(g %in% S), logical(1))
mrk <- function(S) vapply(s5b_ids, function(g){ h <- intersect(g, S)
  if (!length(h)) "" else paste(h, collapse = ",") }, "")
nm  <- function(S) vapply(s5b_ids, function(g){ h <- intersect(g, S)
  if (!length(h)) "" else paste(mo$fuchs_gene_name[match(h, mo$ens)], collapse = ",") }, "")
dr  <- function(S) vapply(s5b_ids, function(g){ h <- intersect(g, S)
  if (!length(h)) "" else { l <- mo$log2FC[match(h[1], mo$ens)]
    sprintf("%s (log2FC=%.2f)", ifelse(l > 0, "up", "down"), l) } }, "")
out <- s5b[, 1:12]
out[["host_gene_ids"]]            <- s5b$gene_ids
out[["n_host_gene_ids"]]          <- s5b$n_gene_ids
out[["MYCN_down_match"]]          <- hit(s5b_ids, DOWN)
out[["MYCN_down_gene_id"]]        <- mrk(DOWN)
out[["MYCN_down_fuchs_name"]]     <- nm(DOWN)
out[["MYCN_down_direction"]]      <- dr(DOWN)
out[["MYCN_up_match_not_tested"]] <- hit(s5b_ids, UPS)
out[["MYCN_up_fuchs_name"]]       <- nm(UPS)
out[["disease"]] <- "Neuroblastoma"
write.csv(out, file.path(OUT, "S5B_36circRNAs_Fuchs2023_MYCN_annotated.csv"), row.names = FALSE)

## ==== overlap table =========================================================
ov <- function(S, set_lbl, tested){
  g <- Reduce(intersect, list(G36, S, BG)); if (!length(g)) return(NULL)
  j <- match(g, mo$ens)
  sy <- vapply(g, function(id) paste(sort(unique(s5b$gene_symbol[
        vapply(s5b_ids, function(v) id %in% v, logical(1))])), collapse = "/"), "")
  data.frame(host_gene = sy, host_gene_id = g, fuchs_gene_name = mo$fuchs_gene_name[j],
             fuchs_set = set_lbl, direction = sprintf("%s (log2FC=%.2f)",
               ifelse(mo$log2FC[j] > 0, "up", "down"), mo$log2FC[j]),
             padj = signif(mo$padj[j], 4), in_fisher_test = tested,
             stringsAsFactors = FALSE) }
ov_gene <- rbind(ov(DOWN, "MNA circs down", TRUE), ov(UPS, "MNA circs up", FALSE))
ov_gene <- ov_gene[order(!ov_gene$in_fisher_test, ov_gene$host_gene), ]
write.csv(ov_gene, file.path(OUT, "overlap_B_hostgene.csv"), row.names = FALSE)

## ==== gene_id / split audit =================================================
aud <- data.frame(circRNA = bsj$circRNA, gene_symbol = bsj$gene_symbol,
                  gene_id_raw = bsj$gene_id,
                  gene_ids_split = vapply(bsj_ids, paste, "", collapse = ","),
                  n_gene_ids = vapply(bsj_ids, length, 1L),
                  in_the36 = bsj$circRNA %in% s5b$circRNA,
                  hits_MYCN_down = vapply(bsj_ids, function(g) any(g %in% DOWN), logical(1)),
                  stringsAsFactors = FALSE)
aud <- aud[order(!aud$in_the36, aud$n_gene_ids < 2, aud$gene_symbol), ]  # the 36 first, then multi-gene
write.csv(aud, file.path(OUT, "id_mapping_audit.csv"), row.names = FALSE)

## ==== UpSet plot ============================================================
Lg <- list(`LOI-36` = intersect(G36, BG), `MYCN-down` = intersect(DOWN, BG))
gp <- grid.grabExpr(print(upset(fromList(Lg), sets = rev(names(Lg)), order.by = "freq",
        keep.order = TRUE, mainbar.y.label = "Shared host genes (Ensembl gene ID)",
        sets.x.label = "host genes per set", text.scale = 1.2)), wrap.grobs = TRUE)
pdf(file.path(OUT, "upset_hostgene.pdf"), width = 8, height = 5)
grid.draw(gp)
invisible(dev.off())

## ==== README ================================================================
multi <- sum(aud$n_gene_ids > 1)
writeLines(c(
"neuroblastoma_Fuchs2023 / MYNC - 36 LOI-circRNA host genes vs MYCN-DOWN circRNAs",
"===============================================================================",
"Script: SCRIPTS/neuroblastoma_Fuchs2023_MYNC.R   (run on rna2 RStudio)",
"Source: Fuchs et al. 2023, Nat Commun, doi:10.1038/s41467-023-38747-4",
"",
"SCOPE",
"  Host-gene level only. The tested category is MOESM6 sheet \"MNA circs down\":",
"  host genes whose circRNA abundance is REDUCED in MYCN-amplified tumours. The",
"  \"MNA circs up\" sheet is read for annotation only and is NOT in the Fisher",
"  test. Coordinate-level matching and the NB-specific category (MOESM8) are out",
"  of scope - see the parent folder.",
"",
"IDENTIFIERS - Ensembl gene ID, no symbol conversion",
"  Fuchs is keyed on its \"ensembl ID\" column (one ID per row, version stripped).",
"  Our side is keyed on the gene_id column of",
"    01_circRNA_DE/limma_BSJ_HnrnpM_vs_NEG_PSD58.tsv",
"  so no gene-symbol -> Ensembl conversion is performed at any point.",
"  (BSJ_FSJ_classification_additive.tsv covers the same 5,820 circRNAs but",
"  carries no gene_id column, which is why the BSJ limma table is read instead.)",
"",
"SPLIT - multi-gene assignments",
"  Some circRNAs are annotated to more than one gene and carry a comma-separated",
"  gene_id, e.g. BBS9 = \"ENSG00000122507,ENSG00000238090\". BOTH the background",
sprintf("  and the 36 are split on those separators (%d of %d circRNAs are affected),", multi, nrow(aud)),
"  so every Ensembl ID is counted once on both sides. Splitting the background",
"  but not the query would shrink the query relative to the background and bias",
"  the odds ratio upward. Fuchs needs no splitting.",
"",
"FISHER UNIVERSE",
"  N = host gene IDs of all HnrnpM-tested circRNAs. Genes are NOT additionally",
"  required to have a circRNA detected in the Fuchs tumour cohort, i.e. MOESM4",
"  is not used and there is no assay-matched/sensitivity pair here.",
"  One- and two-sided Fisher, OR with 95% CI, as in the parent analysis.",
sprintf("  Realised: N=%d  K=%d  n=%d  k=%d  expected=%.2f", N, K, n, k, summ$expected),
sprintf("            OR=%.2f (95%% CI %.2f-%.2f)  P(one-sided)=%.4f  P(two-sided)=%.4f",
        summ$odds_ratio, summ$OR_CI95_low, summ$OR_CI95_high, summ$p_one_sided, summ$p_two_sided),
"  Design check against the agreed target (N=3025 K=296 n=35 k=8 OR=2.78",
sprintf("  CI=1.08-6.36 P=0.0171): %s", if (all(chk)) "PASS" else
        paste("FAIL ->", paste(names(chk)[!chk], collapse = ", "))),
"",
"OUTPUTS",
"  S5B_36circRNAs_Fuchs2023_MYCN_annotated.csv  exact S5B copy (36 rows) + appended:",
"      host_gene_ids, n_host_gene_ids, MYCN_down_match, MYCN_down_gene_id,",
"      MYCN_down_fuchs_name, MYCN_down_direction, MYCN_up_match_not_tested,",
"      MYCN_up_fuchs_name, disease",
"  fisher_summary.csv        the single host-gene test (down-regulated, split)",
"  overlap_B_hostgene.csv    matching host genes; in_fisher_test flags the",
"                            down-regulated ones actually tested",
"  id_mapping_audit.csv      every circRNA, its raw and split gene_id and how many",
"                            Ensembl IDs it carries (the 36 and multi-gene ones first)",
"  upset_hostgene.pdf        UpSet: LOI-36 / MYCN-down (Ensembl gene IDs)",
"  README.txt                this file (written by the script)"),
  file.path(OUT, "README.txt"))

cat(sprintf("MYNC done [down, both split]. N=%d K=%d n=%d k=%d OR=%.2f CI=%.2f-%.2f P=%.4f -- %s\n",
            N, K, n, k, summ$odds_ratio, summ$OR_CI95_low, summ$OR_CI95_high,
            summ$p_one_sided, if (all(chk)) "matches target" else "DOES NOT MATCH TARGET"))
cat(sprintf("  multi-gene circRNAs split: %d of %d ; up-regulated matches (not tested): %d\n",
            multi, nrow(aud), length(Reduce(intersect, list(G36, UPS, BG)))))
