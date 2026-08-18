suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2)
})

# FinalAnalysis item 02 enrichment: GO BP/MF/CC + KEGG ORA on the ADDITIVE gene
# DE list (publication primary in FinalAnalysis). UP, DOWN, ALL stratified.

ROOT     <- "/Volumes/share/mnt/Scratch/PROJECTS/JW_Katana/circRBP_pilot/IV/ForPublication/FinalAnalysis/RESULTS/02_GeneTx_DE"
LIMMA_F  <- file.path(ROOT, "limma_gene_HnrnpM_vs_NEG_PSD58_PSDadd.tsv")
OUT_DIR  <- file.path(ROOT, "Enrichment")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_CUT <- 0.05
LFC_CUT  <- log2(1.5)
TOP_N    <- 20

lm <- read.table(LIMMA_F, sep = "\t", header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE)
lm$ENSEMBL <- sub("\\..*$", "", lm$gene_id)
cat("Universe (additive limma-tested):", nrow(lm), "\n")

sig_mask <- !is.na(lm$adj.P.Val) & lm$adj.P.Val < PADJ_CUT &
            !is.na(lm$logFC) & abs(lm$logFC) >= LFC_CUT
set_UP   <- lm$ENSEMBL[sig_mask & lm$logFC > 0]
set_DOWN <- lm$ENSEMBL[sig_mask & lm$logFC < 0]
set_ALL  <- lm$ENSEMBL[sig_mask]
universe <- lm$ENSEMBL
cat(sprintf("UP: %d, DOWN: %d, ALL: %d (adj.P<%.2f & |FC|>=%.2f)\n",
            length(set_UP), length(set_DOWN), length(set_ALL), PADJ_CUT, 1.5))

suppressMessages({
  e2e_df <- bitr(universe, fromType = "ENSEMBL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = FALSE)
})
e2e_df <- e2e_df[!is.na(e2e_df$ENTREZID), ]
e2e_df <- e2e_df[!duplicated(e2e_df$ENSEMBL), ]
e2e <- setNames(e2e_df$ENTREZID, e2e_df$ENSEMBL)
ent_UP   <- unique(e2e[set_UP][!is.na(e2e[set_UP])])
ent_DOWN <- unique(e2e[set_DOWN][!is.na(e2e[set_DOWN])])
ent_ALL  <- unique(e2e[set_ALL][!is.na(e2e[set_ALL])])
ent_universe <- unique(e2e[!is.na(e2e)])

run_go <- function(genes, set_lab, ont) {
  if (length(genes) < 5) return(NULL)
  res <- tryCatch(
    enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
             keyType = "ENSEMBL", ont = ont, pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
    error = function(e) NULL)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  df <- as.data.frame(res)
  write.table(df, file.path(OUT_DIR, sprintf("HnrnpM_%s_GO_%s.tsv", set_lab, ont)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(file.path(OUT_DIR, sprintf("HnrnpM_%s_GO_%s.pdf", set_lab, ont)),
         dotplot(res, showCategory = TOP_N) +
           ggtitle(sprintf("HnrnpM %s - GO %s (additive)", set_lab, ont)),
         width = 9, height = 7)
  list(df = df, res = res)
}

run_kegg <- function(genes, ent_universe, set_lab) {
  if (length(genes) < 5) return(NULL)
  res <- tryCatch(
    enrichKEGG(gene = genes, universe = ent_universe, organism = "hsa",
               keyType = "kegg", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) NULL)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  df <- as.data.frame(res)
  write.table(df, file.path(OUT_DIR, sprintf("HnrnpM_%s_KEGG.tsv", set_lab)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(file.path(OUT_DIR, sprintf("HnrnpM_%s_KEGG.pdf", set_lab)),
         dotplot(res, showCategory = TOP_N) +
           ggtitle(sprintf("HnrnpM %s - KEGG (additive)", set_lab)),
         width = 9, height = 7)
  list(df = df, res = res)
}

sets_e <- list(UP = set_UP, DOWN = set_DOWN, ALL = set_ALL)
sets_n <- list(UP = ent_UP, DOWN = ent_DOWN, ALL = ent_ALL)
summary_rows <- list()
for (s in names(sets_e)) {
  cat("---", s, "---\n")
  for (ont in c("BP","MF","CC")) {
    r <- run_go(sets_e[[s]], s, ont)
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      set = s, source = paste0("GO_", ont),
      n_input = length(sets_e[[s]]),
      n_sig_terms = if (is.null(r)) 0L else nrow(r$df),
      top_term = if (is.null(r)) NA else r$df$Description[1],
      top_padj = if (is.null(r)) NA else r$df$p.adjust[1])
  }
  r <- run_kegg(sets_n[[s]], ent_universe, s)
  summary_rows[[length(summary_rows)+1]] <- data.frame(
    set = s, source = "KEGG",
    n_input = length(sets_n[[s]]),
    n_sig_terms = if (is.null(r)) 0L else nrow(r$df),
    top_term = if (is.null(r)) NA else r$df$Description[1],
    top_padj = if (is.null(r)) NA else r$df$p.adjust[1])
}
summ <- do.call(rbind, summary_rows)
write.table(summ, file.path(OUT_DIR, "enrichment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)
cat("\nDone.\n")
