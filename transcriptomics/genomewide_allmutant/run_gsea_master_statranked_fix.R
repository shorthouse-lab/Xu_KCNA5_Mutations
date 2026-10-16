## ---------------------------------------------------------------------------
## PROJECT-WIDE FIX: every GSEA table in this project ranked genes by
## sign(log2FoldChange) * -log10(padj) -- a statistic that collapses ~70% of
## genes into tied values given this cohort's padj distribution (only ~5,900
## unique padj values across ~28,000 tested genes). fgsea's own documentation
## warns tied preranked stats make results order-dependent/arbitrary; this
## was verified to be a REAL problem, not boilerplate caution -- independent
## reruns of the identical code on identical inputs produced different NES,
## including full sign flips for some pathways (e.g. REACTOME_G2_M_
## CHECKPOINTS and the entire AllMutant_vs_WT Hallmark table). Discovered
## while fact-checking a manuscript paragraph's exact pathway statistics.
##
## Fix: rank by DESeq2's own Wald `stat` column instead -- continuous, zero
## ties (one unique value per gene). Verified stable across 6 different
## random seeds after the fix (previously flipped sign between runs).
##
## Regenerates every fgsea_<contrast>_<collection>.csv this project uses,
## for all 4 collections (Hallmark/KEGG/Reactome/Oncogenic), across every
## contrast with an existing GSEA table:
##   - AllMutant_vs_WT
##   - KCNA5_Low_vs_High_WT
##   - T1_vs_WT, S1S4_vs_WT (domain-split)
##   - T1_vs_S1S4 (direct domain contrast)
## Overwrites the same canonical filenames every downstream script already
## reads, so no other script needs path changes -- only re-running.
## ---------------------------------------------------------------------------
library(fgsea)
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS_MAIN <- "analysis3_transcriptomic"
ANALYSIS_DOMAIN <- "analysis3_domain_split"

build_ranks_stat <- function(df) {
  df <- df %>% filter(!is.na(stat), !is.na(gene_symbol), gene_symbol != "")
  df <- df %>% group_by(gene_symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  ranks <- setNames(df$stat, df$gene_symbol)
  sort(ranks, decreasing = TRUE)
}

run_fgsea_ml <- function(ranks, pathways, label) {
  log_progress("  fgseaMultilevel (", label, "): ", length(ranks), " genes (", length(unique(ranks)),
               " unique), ", length(pathways), " pathways...")
  set.seed(42)
  res <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500, eps = 0)
  res$contrast <- label
  res %>% arrange(padj)
}

save_fn <- function(res, name, analysis) {
  out <- res %>% mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ","))
  save_table(as.data.frame(out), name, analysis)
}

collections <- list(
  hallmark  = "hallmark_genesets.rds",
  kegg      = "kegg_genesets.rds",
  reactome  = "reactome_genesets.rds",
  oncogenic = "oncogenic_genesets.rds"
)
pathway_sets <- lapply(collections, function(f) readRDS(file.path(RAW_DIR, f)))

## =============================================================================
## A. AllMutant_vs_WT (analysis3_transcriptomic)
## =============================================================================
log_progress("=== AllMutant_vs_WT ===")
de_allmut <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS_MAIN, "deseq2_AllMutant_vs_WT.csv"))
ranks_allmut <- build_ranks_stat(de_allmut)
for (coll in names(collections)) {
  res <- run_fgsea_ml(ranks_allmut, pathway_sets[[coll]], "AllMutant_vs_WT")
  log_progress("  -> ", coll, ": ", sum(res$padj < 0.25, na.rm = TRUE), " at FDR q<0.25, ",
               sum(res$padj < 0.05, na.rm = TRUE), " at q<0.05")
  save_fn(res, paste0("fgsea_AllMutant_vs_WT_", coll), ANALYSIS_MAIN)
}

## =============================================================================
## B. KCNA5_Low_vs_High_WT (analysis3_transcriptomic)
## =============================================================================
log_progress("=== KCNA5_Low_vs_High_WT ===")
de_expr <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS_MAIN, "deseq2_KCNA5_Low_vs_High_WT.csv"))
ranks_expr <- build_ranks_stat(de_expr)
for (coll in names(collections)) {
  res <- run_fgsea_ml(ranks_expr, pathway_sets[[coll]], "KCNA5_Low_vs_High_WT")
  log_progress("  -> ", coll, ": ", sum(res$padj < 0.25, na.rm = TRUE), " at FDR q<0.25")
  save_fn(res, paste0("fgsea_KCNA5_Low_vs_High_WT_", coll), ANALYSIS_MAIN)
}

## =============================================================================
## C. Domain-split: T1_vs_WT, S1S4_vs_WT, T1_vs_S1S4 (analysis3_domain_split)
## =============================================================================
log_progress("=== Domain-split (T1_vs_WT, S1S4_vs_WT, T1_vs_S1S4) ===")
res_list <- readRDS(file.path(PROCESSED_DIR, "analysis3_domain_deseq2_results.rds"))
ranks_T1   <- build_ranks_stat(res_list$res_T1_vs_WT)
ranks_S1S4 <- build_ranks_stat(res_list$res_S1S4_vs_WT)
ranks_T1S4 <- build_ranks_stat(res_list$res_T1_vs_S1S4)

for (coll in names(collections)) {
  fg_T1   <- run_fgsea_ml(ranks_T1,   pathway_sets[[coll]], "T1_vs_WT")
  fg_S1S4 <- run_fgsea_ml(ranks_S1S4, pathway_sets[[coll]], "S1S4_vs_WT")
  fg_T1S4 <- run_fgsea_ml(ranks_T1S4, pathway_sets[[coll]], "T1_vs_S1S4")

  if (coll == "hallmark") {
    save_fn(fg_T1, "fgsea_T1_vs_WT_domain", ANALYSIS_DOMAIN)
    save_fn(fg_S1S4, "fgsea_S1S4_vs_WT_domain", ANALYSIS_DOMAIN)
    save_fn(fg_T1S4, "fgsea_T1_vs_S1S4_domain", ANALYSIS_DOMAIN)
  } else {
    save_fn(fg_T1, paste0("fgsea_T1_vs_WT_", coll), ANALYSIS_DOMAIN)
    save_fn(fg_S1S4, paste0("fgsea_S1S4_vs_WT_", coll), ANALYSIS_DOMAIN)
    save_fn(fg_T1S4, paste0("fgsea_T1_vs_S1S4_", coll), ANALYSIS_DOMAIN)
  }
  log_progress("  -> ", coll, ": T1 q<0.25=", sum(fg_T1$padj < 0.25, na.rm = TRUE),
               ", S1S4 q<0.25=", sum(fg_S1S4$padj < 0.25, na.rm = TRUE),
               ", T1vsS1S4 q<0.25=", sum(fg_T1S4$padj < 0.25, na.rm = TRUE))

  merged <- inner_join(
    fg_T1 %>% transmute(pathway, pval_T1 = pval, NES_T1 = NES, padj_T1 = padj),
    fg_S1S4 %>% transmute(pathway, pval_S1S4 = pval, NES_S1S4 = NES, padj_S1S4 = padj),
    by = "pathway"
  )
  save_table(merged, paste0(coll, "_T1_S1S4_merged_for_barplot"), ANALYSIS_DOMAIN)
}

log_progress("Master stat-ranked GSEA regeneration complete.")
