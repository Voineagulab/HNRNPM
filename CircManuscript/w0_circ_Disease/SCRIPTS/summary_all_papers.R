#!/usr/bin/env Rscript
# =============================================================================
# summary_all_papers.R
# Consolidates the four disease-overlap analyses (Fuchs/NB, Song/glioma,
# Dube/AD, Phillips/AD) into ONE workbook:
#   RESULTS/circDE36_disease_overlap_summary.xlsx
# Sheets:
#   * one per paper (the annotated S5B tables)
#   * circRNA_overlap    (one row per paper, like the clean workbook)
#   * Notes
#   * per_circRNA_disease (36 rows x disease; "Implicated in ...? (match level)"
#     columns highlighted GREEN=enriched / RED=depleted; a separate
#     "...: enriched/depleted?" text column; blank/— for detected-not-significant)
# Zhao 2024 is excluded (juncid not resolvable to coordinates/genes); the three
# psychiatric/schizophrenia papers are excluded per project scope.
# Run on rna2 RStudio.
# =============================================================================
suppressPackageStartupMessages({ library(openxlsx) })
RES <- "/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/JW/CircManuscript/w0_circ_Disease/RESULTS"
rd  <- function(p) read.csv(file.path(RES, p), check.names = FALSE, stringsAsFactors = FALSE)
tf  <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE","T","1")
fu <- rd("neuroblastoma_Fuchs2023/S5B_36circRNAs_Fuchs2023_annotated.csv")
so <- rd("glioblastoma_Song2016/S5B_36circRNAs_Song2016_annotated.csv")
du <- rd("AD_Dube2019/S5B_36circRNAs_Dube2019_annotated.csv")
ph <- rd("AD_Phillips2026/S5B_36circRNAs_Phillips2026_annotated.csv")

# canonical 36 order/keys from Fuchs table (all share circRNA + gene_symbol)
circ <- fu$circRNA; gene <- fu$gene_symbol; n36 <- length(circ)
idx <- function(df) match(circ, df$circRNA)
dir_word <- function(s){ s <- tolower(as.character(s))
  ifelse(grepl("up", s) & !grepl("down", s), "Enriched",
  ifelse(grepl("down", s), "Depleted", "")) }

## ---- per-disease implication (match level + direction) ---------------------
# Alzheimer = Dube (host-gene) + Phillips (coordinate/host-gene)
iph <- idx(ph); idu <- idx(du)
ph_A <- tf(ph$AD_A_coordinate_match)[iph]; ph_B <- tf(ph$AD_B_hostgene_match)[iph]
du_B <- tf(du$AD_B_hostgene_match)[idu]
AD_impl  <- ph_A | ph_B | du_B
AD_level <- ifelse(ph_A, "COORDINATE", ifelse(ph_B | du_B, "host-gene", ""))
AD_dirsrc<- ifelse(ph_A, ph$AD_A_direction[iph], ifelse(ph_B, ph$AD_B_direction[iph],
             ifelse(du_B, du$AD_B_direction[idu], "")))
AD_dir   <- dir_word(AD_dirsrc)

# Neuroblastoma = Fuchs (NB-specific coord/gene, prioritised over MYCN host-gene)
f_A <- tf(fu$NBspecific_A_coordinate_match); f_Bs <- tf(fu$NBspecific_B_hostgene_match)
f_Bm <- tf(fu$MYCN_B_hostgene_match)
NB_impl  <- f_A | f_Bs | f_Bm
NB_level <- ifelse(f_A, "COORDINATE", ifelse(f_Bs | f_Bm, "host-gene", ""))
NB_dirsrc<- ifelse(f_A | f_Bs, "up (NB-specific)", ifelse(f_Bm, fu$MYCN_B_direction, ""))
NB_dir   <- dir_word(NB_dirsrc)

# Glioblastoma = Song (coordinate q<0.05 / host-gene)
iso <- idx(so); g_A <- tf(so$glioma_A_coordinate_match)[iso]; g_B <- tf(so$glioma_B_hostgene_match)[iso]
GB_impl  <- g_A | g_B
GB_level <- ifelse(g_A, "COORDINATE", ifelse(g_B, "host-gene", ""))
GB_dirsrc<- ifelse(g_A, so$glioma_A_direction[iso], ifelse(g_B, so$glioma_B_direction[iso], ""))
GB_dir   <- dir_word(GB_dirsrc)

cellval <- function(impl, level) ifelse(impl, paste0("YES - ", level), "—")
pcd <- data.frame(
  `circRNA (hg38 coordinate)` = circ, `Host gene` = gene,
  `Implicated in Alzheimer's? (match level)` = cellval(AD_impl, AD_level),
  `Alzheimer: enriched/depleted?` = AD_dir,
  `Implicated in neuroblastoma? (match level)` = cellval(NB_impl, NB_level),
  `Neuroblastoma: enriched/depleted?` = NB_dir,
  `Implicated in glioblastoma? (match level)` = cellval(GB_impl, GB_level),
  `Glioblastoma: enriched/depleted?` = GB_dir,
  check.names = FALSE, stringsAsFactors = FALSE)

## ---- circRNA_overlap (one row per paper) -----------------------------------
overlap <- data.frame(
  `Research title` = c(
    "An atlas of cortical circular RNA expression in Alzheimer disease brains",
    "Blood-based circular RNAs for early diagnosis of Alzheimer's disease",
    "Defining the landscape of circular RNAs in neuroblastoma (MYCN)",
    "Circular RNA profile in gliomas revealed by identification tool UROBORUS"),
  `Citation` = c("(Dube et al., 2019)","(Phillips et al., 2026)","(Fuchs et al., 2023)","(Song et al., 2016)"),
  `Disease` = c("Alzheimer","Alzheimer","Neuroblastoma","Glioblastoma"),
  `Overlapping circRNAs (level specified)` = c(
    "HOST-GENE only (no coordinates in Dube tables): circATRNL1, circPTK2 (2/32 host genes; AD_any).",
    "COORDINATE + host-gene: circSCLT1, circRBM33 (AD-associated FDR<0.05, up in AD blood).",
    "COORDINATE (NB-specific, MOESM8): circEML5, circZNF800 (up). HOST-GENE (MYCN-DE, MOESM6): ADAM22, BBS9, FBXL13, KANSL1L, MAP3K4, TBCD, ZNF800.",
    "COORDINATE (q<0.05, recomputed): 7 circRNAs (circEML5, circZBTB44, circZNF800, circATRNL1, circCSNK1G3, circMAP3K4, circTNPO3) all DEPLETED in GBM; none enriched."),
  `Fisher (primary universe)` = c(
    "host-gene OR=1.38, p=0.44 (Knight_3547)",
    "coordinate OR=2.25, p=0.27 (assay-matched)",
    "coordinate OR=11.56, p=0.017 (assay-matched); MYCN host-gene OR=1.83, p=0.14",
    "coordinate OR=0.33, p=0.98 (NOT enriched; glioma has global circRNA loss)"),
  `Supplementary table(s)` = c("SuppTable 11/13/15","Table S3/S4 (+S5)","MOESM8/MOESM6 (+MOESM4 universe)","Table S5 (+S4 groups; DE recomputed)"),
  `Script / folder` = c("RESULTS/AD_Dube2019.R","RESULTS/AD_Phillips2026.R","RESULTS/neuroblastoma_Fuchs2023.R","RESULTS/glioblastoma_Song2016.R"),
  check.names = FALSE, stringsAsFactors = FALSE)

notes <- data.frame(Notes = c(
  "Summary of the overlap between the 36 HNRNPM-regulated (independent-backsplicing) circRNAs and disease-associated circRNAs from four public datasets.",
  "Reference set: sheet S5B_independent_circRNAs of Supplementary_Tables_JW_2.xlsx (36 circRNAs, hg38). Universe: BSJ_FSJ_classification_additive.tsv (5,820 HnrnpM-tested circRNAs).",
  "(A) coordinate match = exact hg38 back-splice junction, +/-1 nt. (B) host-gene match = gene symbol.",
  "Disease-associated = statistically significant in each paper: Phillips FDR<0.05 (S3/S4); Song q<0.05 (recomputed Wilcoxon rank-sum GBM vs normal + BH from Table S5/S4); Fuchs NB-specific (MOESM8) + MYCN-amplified DE (MOESM6); Dube listed (all FDR<0.05).",
  "Song/glioma: the paper publishes no fold change; disease-association and direction are OUR recomputation from Supplementary Table S5 (see RESULTS/glioblastoma_Song2016/README.txt). Fold changes shown for glioma are DERIVED, not published.",
  "'Merely detected but not significant' circRNAs are used only to build the Fisher universe and are NOT flagged A/B; they appear blank/dash in per_circRNA_disease.",
  "per_circRNA_disease: GREEN = enriched (up in disease), RED = depleted (down in disease), in the 'Implicated in ...? (match level)' columns; dash = not implicated.",
  "EXCLUSIONS: Zhao et al. 2024 (circMeta2) dropped - its circRNAs are identified only by an opaque integer 'juncid' with no coordinate or gene mapping in the deposited tables. The three schizophrenia/psychiatric papers are out of scope for this summary.",
  "Fisher tests are one-sided (greater) for enrichment; full statistics per paper in each RESULTS/<paper>/fisher_summary.csv."),
  check.names = FALSE, stringsAsFactors = FALSE)

## ---- write workbook with styling -------------------------------------------
wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#1F4E5F", fontColour = "white",
                   halign = "left", valign = "top", wrapText = TRUE, border = "TopBottomLeftRight")
wrapS <- createStyle(wrapText = TRUE, valign = "top")
greenS <- createStyle(fgFill = "#C6EFCE", fontColour = "#1B7A43", textDecoration = "bold", wrapText = TRUE, valign = "top")
redS   <- createStyle(fgFill = "#FFC7CE", fontColour = "#B24B3A", textDecoration = "bold", wrapText = TRUE, valign = "top")
addsheet <- function(name, df, widths){ addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr); freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_along(df), widths = widths)
  addStyle(wb, name, wrapS, rows = 2:(nrow(df)+1), cols = seq_along(df), gridExpand = TRUE, stack = TRUE) }

# per-paper sheets (annotated tables)
addsheet("Fuchs2023_NB", fu, c(26,10, rep(9,10), rep(16,9))[seq_len(ncol(fu))])
addsheet("Song2016_glioma", so, c(26,10, rep(9,10), rep(20,5))[seq_len(ncol(so))])
addsheet("Dube2019_AD", du, c(26,10, rep(9,10), rep(20,5))[seq_len(ncol(du))])
addsheet("Phillips2026_AD", ph, c(26,10, rep(9,10), rep(22,5))[seq_len(ncol(ph))])

# circRNA_overlap + Notes
addsheet("circRNA_overlap", overlap, c(42,20,14,60,42,26,26))
addWorksheet(wb, "Notes"); writeData(wb, "Notes", notes, headerStyle = hdr)
setColWidths(wb, "Notes", cols = 1, widths = 120)
addStyle(wb, "Notes", wrapS, rows = 2:(nrow(notes)+1), cols = 1, gridExpand = TRUE, stack = TRUE)

# per_circRNA_disease with green/red on the 3 match-level columns
addsheet("per_circRNA_disease", pcd, c(27,12,31,20,31,20,31,20))
lvl_cols <- c(3,5,7); dir_vec <- list(AD_dir, NULL, NB_dir, NULL, GB_dir)
for (j in lvl_cols){ dcol <- switch(as.character(j), "3"=AD_dir, "5"=NB_dir, "7"=GB_dir)
  for (i in seq_len(n36)){
    if (dcol[i] == "Enriched") addStyle(wb, "per_circRNA_disease", greenS, rows = i+1, cols = j, stack = TRUE)
    else if (dcol[i] == "Depleted") addStyle(wb, "per_circRNA_disease", redS, rows = i+1, cols = j, stack = TRUE) } }

out <- file.path(RES, "circDE36_disease_overlap_summary.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
cat(sprintf("Wrote %s\n  AD implicated=%d ; NB implicated=%d ; GBM implicated=%d\n",
            out, sum(AD_impl), sum(NB_impl), sum(GB_impl)))
