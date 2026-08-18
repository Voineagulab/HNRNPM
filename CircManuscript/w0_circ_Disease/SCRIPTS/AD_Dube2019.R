#!/usr/bin/env Rscript
# =============================================================================
# AD_Dube2019.R
# Overlap of the 36 HNRNPM-regulated circRNAs with Alzheimer's-disease-associated
# circRNAs from Dube et al. 2019 (Nat Neurosci).
# HOST-GENE LEVEL ONLY: Dube collapse back-splice counts onto the host gene and
# report "circGENE" + chromosome only (no start/end), so coordinate matching (A)
# is not possible (column A = FALSE for all; A_direction = NA). Column B = host
# gene symbol match.
# AD-associated sets = SuppTable 11 (CDR meta), 13 (Braak meta), 15 (case-control
# meta); every circRNA listed passed the paper's FDR<0.05. Direction from log2FC.
# Outputs -> RESULTS/AD_Dube2019/
# Run on rna2 RStudio.
# =============================================================================
suppressPackageStartupMessages({ library(readxl); library(UpSetR) })

BASE  <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript"
W     <- file.path(BASE, "w0_circ_Disease")
S5BX  <- file.path(BASE, "Supplementary_Table", "Supplementary_Tables_JW_2.xlsx")
IVCLS <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/01_circRNA_DE/FSJ/BSJ_FSJ_classification_additive.tsv"
XLS   <- list.files(file.path(W, "Alzheimer", "Dube2019"), pattern = "sup.xls$", full.names = TRUE)[1]
OUT   <- file.path(W, "RESULTS", "AD_Dube2019")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

up <- function(x) toupper(trimws(as.character(x)))
clean_gene <- function(x) toupper(sub("^circ", "", trimws(as.character(x))))
frow <- function(category, level, universe, k, n, K, N, note = "") {
  if (is.na(N) || N <= 0 || (N-K-(n-k)) < 0) return(data.frame(category, match_level = level, universe, N = N,
    K_category = K, n_HnrnpM = n, k_overlap = k, expected = NA, odds_ratio = NA, OR_CI95_low = NA,
    OR_CI95_high = NA, p_one_sided = NA, p_two_sided = NA, note = paste(note, "invalid table"), stringsAsFactors = FALSE))
  tab <- matrix(c(k, K-k, n-k, N-K-(n-k)), 2, byrow = TRUE)
  f1 <- fisher.test(tab, alternative = "greater"); f2 <- fisher.test(tab)
  data.frame(category, match_level = level, universe, N = N, K_category = K, n_HnrnpM = n, k_overlap = k,
    expected = round(n*K/N, 3), odds_ratio = round(unname(f2$estimate), 2),
    OR_CI95_low = round(f2$conf.int[1], 2), OR_CI95_high = round(f2$conf.int[2], 2),
    p_one_sided = signif(f1$p.value, 4), p_two_sided = signif(f2$p.value, 4), note, stringsAsFactors = FALSE)
}

## ---- the 36 (host genes) + HnrnpM universe genes ---------------------------
s5b <- as.data.frame(read_excel(S5BX, sheet = "S5B_independent_circRNAs", skip = 1))
s5b <- s5b[!is.na(s5b$circRNA) & s5b$circRNA != "", ]; s5b$gene <- up(s5b$gene_symbol)
cls <- read.delim(IVCLS, stringsAsFactors = FALSE, check.names = FALSE); cls$gene <- up(cls$gene_symbol)
g36  <- sort(unique(s5b$gene[s5b$gene != ""]))
hn_genes <- sort(unique(cls$gene[cls$gene != ""]))

## ---- Dube meta tables (host-gene + log2FC direction) -----------------------
read_meta <- function(sheet){
  raw <- suppressMessages(as.data.frame(read_excel(XLS, sheet = sheet, col_names = FALSE)))
  hr  <- which(apply(raw, 1, function(r) any(trimws(as.character(r)) == "circRNA")))[1]
  hdr <- trimws(as.character(unlist(raw[hr, ])))
  c_gene <- which(hdr == "circRNA")[1]; c_disc <- which(hdr == "log2FC")[1]
  dat <- raw[(hr+1):nrow(raw), , drop = FALSE]; circ <- trimws(as.character(dat[[c_gene]]))
  keep <- !is.na(circ) & grepl("^circ", circ); dat <- dat[keep, , drop = FALSE]; circ <- circ[keep]
  data.frame(circRNA = circ, gene = clean_gene(circ),
             disc_log2FC = suppressWarnings(as.numeric(dat[[c_disc]])), stringsAsFactors = FALSE)
}
cdr <- read_meta("SuppTable11"); braak <- read_meta("SuppTable13"); cc <- read_meta("SuppTable15")
AD_sets <- lapply(list(CDR = cdr$gene, Braak = braak$gene, case_control = cc$gene),
                  function(g) sort(unique(g[g != ""])))
AD_sets$AD_any <- sort(unique(unlist(AD_sets)))
# per-gene Dube direction: prefer case-control DE, else CDR, else Braak
dube_l2 <- function(g){ v <- cc$disc_log2FC[match(g, cc$gene)]
  v[is.na(v)] <- cdr$disc_log2FC[match(g[is.na(v)], cdr$gene)]
  v[is.na(v)] <- braak$disc_log2FC[match(g[is.na(v)], braak$gene)]; v }
# which AD traits each gene hits (for the direction annotation)
traits <- function(g) vapply(g, function(x) paste(c("CDR","Braak","case_control")[c(
  x %in% AD_sets$CDR, x %in% AD_sets$Braak, x %in% AD_sets$case_control)], collapse = "+"), "")

## ---- universe of all Dube-reported host genes (union of all supp tables) ----
all_dube_genes <- local({ sh <- excel_sheets(XLS); s <- character(0)
  for (x in sh){ raw <- tryCatch(suppressMessages(as.data.frame(read_excel(XLS, sheet = x, col_names = FALSE))),
                                 error = function(e) NULL)
    if (is.null(raw) || !ncol(raw)) next; col1 <- trimws(as.character(raw[[1]]))
    s <- c(s, clean_gene(col1[grepl("^circ", col1)])) }; sort(unique(s)) })
NCOUNT <- c(Knight_discovery_3547 = 3547, MSBB_BM44_4330 = 4330, all_human_circ_12000 = 12000)

## ==== S5B CSV with appended A/B + direction (A not assessable) ==============
B <- s5b$gene %in% AD_sets$AD_any
Bdir <- ifelse(B, sprintf("%s in AD (log2FC=%.2f; %s)", ifelse(dube_l2(s5b$gene) > 0, "up", "down"),
                          dube_l2(s5b$gene), traits(s5b$gene)), "")
out <- s5b[, 1:12]
out[["AD_A_coordinate_match"]] <- rep(FALSE, nrow(s5b))
out[["AD_A_direction"]]        <- rep("NA (Dube host-gene + chromosome only, no coordinates)", nrow(s5b))
out[["AD_B_hostgene_match"]]   <- B
out[["AD_B_direction"]]        <- Bdir
out[["disease"]] <- "Alzheimer"
write.csv(out, file.path(OUT, "S5B_36circRNAs_Dube2019_annotated.csv"), row.names = FALSE)

## ==== Fisher (host-gene; AD_any + per-trait; multiple universes) ============
n <- length(g36)
rows <- list()
for (nm in names(NCOUNT)) rows[[length(rows)+1]] <- frow("AD-associated (AD_any)", "host-gene",
  sprintf("count-based: %s (N=%d)", nm, NCOUNT[[nm]]),
  length(intersect(g36, AD_sets$AD_any)), n, length(AD_sets$AD_any), NCOUNT[[nm]])
# transparency-only reported-genes universe
rows[[length(rows)+1]] <- frow("AD-associated (AD_any)", "host-gene",
  sprintf("Dube_reported_only_biased (N=%d) [NOT for inference]", length(all_dube_genes)),
  length(Reduce(intersect, list(g36, AD_sets$AD_any, all_dube_genes))),
  length(intersect(g36, all_dube_genes)), length(intersect(AD_sets$AD_any, all_dube_genes)),
  length(all_dube_genes))
for (s in c("CDR","Braak","case_control")) rows[[length(rows)+1]] <- frow(
  sprintf("AD-associated (%s)", s), "host-gene", "count-based: Knight_discovery_3547 (N=3547) [PRIMARY]",
  length(intersect(g36, AD_sets[[s]])), n, length(AD_sets[[s]]), 3547)
summ <- do.call(rbind, rows)
write.csv(summ, file.path(OUT, "fisher_summary.csv"), row.names = FALSE)

## ==== overlap table =========================================================
gi <- intersect(g36, AD_sets$AD_any)
write.csv(data.frame(host_gene = gi, AD_direction = ifelse(dube_l2(gi) > 0, "up in AD", "down in AD"),
                     log2FC = round(dube_l2(gi), 3), traits = traits(gi), stringsAsFactors = FALSE),
          file.path(OUT, "overlap_B_hostgene.csv"), row.names = FALSE)

## ==== UpSet (host-gene, PDF) ================================================
L <- list(`HnrnpM-36` = g36, CDR = AD_sets$CDR, Braak = AD_sets$Braak, `case-control` = AD_sets$case_control)
pdf(file.path(OUT, "upset_hostgene.pdf"), width = 8, height = 5)
print(upset(fromList(L), sets = rev(names(L)), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Shared host genes", sets.x.label = "genes per set", text.scale = 1.2))
dev.off()

cat(sprintf("Dube done. AD_any genes=%d; overlap with 36 host genes=%d (%s)\n",
            length(AD_sets$AD_any), length(gi), paste(gi, collapse = ", ")))
