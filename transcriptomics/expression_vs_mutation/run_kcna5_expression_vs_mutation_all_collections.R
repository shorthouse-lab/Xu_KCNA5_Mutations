## ---------------------------------------------------------------------------
## Extends run_deseq2_kcna5_expression_groups.R's Hallmark-only concordance
## check (KCNA5-low-vs-high in WT patients vs AllMutant-vs-WT) to the other
## three pathway collections (KEGG, Reactome, oncogenic signatures), so the
## "does mutation phenocopy low expression?" question is answered across all
## pathway collections used elsewhere in this study, not just Hallmark.
##
## GSEA computation itself lives in run_gsea_master_statranked_fix.R
## (tie-free ranking by DESeq2's Wald `stat` column, covering this contrast
## for all 4 collections) -- run that script first; it writes
## fgsea_KCNA5_Low_vs_High_WT_{hallmark,kegg,reactome,oncogenic}.csv. This
## script only loads those already-computed tables and builds the
## concordance comparison against AllMutant_vs_WT -- no fgsea call lives
## here, so there is no ranking logic to go stale.
## ---------------------------------------------------------------------------
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_transcriptomic"

collections <- list(
  kegg = list(expr_file = "fgsea_KCNA5_Low_vs_High_WT_kegg.csv", allmut_file = "fgsea_AllMutant_vs_WT_kegg.csv"),
  reactome = list(expr_file = "fgsea_KCNA5_Low_vs_High_WT_reactome.csv", allmut_file = "fgsea_AllMutant_vs_WT_reactome.csv"),
  oncogenic = list(expr_file = "fgsea_KCNA5_Low_vs_High_WT_oncogenic.csv", allmut_file = "fgsea_AllMutant_vs_WT_oncogenic.csv")
)

all_concordance <- list()
for (coll_name in names(collections)) {
  spec <- collections[[coll_name]]
  expr_path <- file.path(RESULTS_DIR, "tables", ANALYSIS, spec$expr_file)
  if (!file.exists(expr_path)) {
    stop("Missing ", expr_path, " -- run run_gsea_master_statranked_fix.R first.")
  }
  log_progress("--- ", coll_name, ": loading pre-computed KCNA5_Low_vs_High_WT fgsea ---")
  fg_expr_df <- read.csv(expr_path)

  allmut_fg <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, spec$allmut_file))
  sig_mut_pathways <- allmut_fg %>% filter(padj < 0.05) %>% pull(pathway)

  pathway_merge <- inner_join(
    fg_expr_df %>% filter(pathway %in% sig_mut_pathways) %>% transmute(pathway, NES_expr = NES, padj_expr = padj),
    allmut_fg %>% filter(pathway %in% sig_mut_pathways) %>% transmute(pathway, NES_mut = NES, padj_mut = padj),
    by = "pathway"
  ) %>% mutate(collection = coll_name, concordant = sign(NES_expr) == sign(NES_mut)) %>% arrange(NES_mut)

  save_table(pathway_merge, paste0("kcna5_expression_vs_mutation_", coll_name, "_concordance"), ANALYSIS)
  log_progress(coll_name, ": ", length(sig_mut_pathways), " pathways significant in AllMutant_vs_WT, ",
               sum(pathway_merge$concordant), "/", nrow(pathway_merge), " concordant in direction with KCNA5-low-vs-high")
  all_concordance[[coll_name]] <- pathway_merge
}

## =============================================================================
## Combined summary across all 4 collections (Hallmark + KEGG + Reactome + oncogenic)
## =============================================================================
hallmark_concordance <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "kcna5_expression_vs_mutation_hallmark_concordance.csv")) %>%
  mutate(collection = "hallmark", concordant = sign(NES_expr) == sign(NES_mut))

combined <- bind_rows(hallmark_concordance, all_concordance$kegg, all_concordance$reactome, all_concordance$oncogenic)
save_table(combined, "kcna5_expression_vs_mutation_all_collections_concordance", ANALYSIS)

summary_by_collection <- combined %>% group_by(collection) %>%
  summarise(n_pathways = n(), n_concordant = sum(concordant), n_discordant = sum(!concordant), .groups = "drop")
log_progress("\n=== Combined concordance summary across all 4 pathway collections ===\n",
             paste(capture.output(print(summary_by_collection)), collapse = "\n"))
log_progress("\nOverall: ", sum(combined$concordant), "/", nrow(combined),
             " significant AllMutant_vs_WT pathways show concordant direction with KCNA5-low-vs-high (WT)")
log_progress("\nDiscordant pathways (opposite direction):\n",
             paste(capture.output(print(combined %>% filter(!concordant) %>% select(collection, pathway, NES_expr, NES_mut), digits=3)), collapse = "\n"))

log_progress("All-collections KCNA5 expression vs mutation concordance analysis complete.")
