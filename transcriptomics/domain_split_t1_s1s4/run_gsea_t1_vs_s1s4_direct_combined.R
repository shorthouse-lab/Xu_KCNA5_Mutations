## ---------------------------------------------------------------------------
## Direct T1-vs-S1S4 GSEA, all 4 collections pooled into one figure.
##
## GSEA computation itself lives entirely in run_gsea_master_statranked_fix.R
## (tie-free ranking by DESeq2's Wald `stat` column -- see that script's
## header for why the older sign(log2FC)*-log10(padj) ranking was replaced).
## Run that script FIRST; it writes fgsea_T1_vs_S1S4_domain.csv (Hallmark)
## and fgsea_T1_vs_S1S4_{kegg,reactome,oncogenic}.csv. This script only pools
## those 4 already-computed tables into one figure -- it does not run fgsea
## itself, so there is no ranking logic here to go stale.
##
## "Excluding those with only a few genes driving": a GSEA pathway can be
## nominally significant while its leading edge -- the actual genes driving
## the enrichment score -- is tiny (2-3 genes), which is exactly the
## thin-evidence problem already flagged for the GO BP ORA result on this
## contrast (57 terms, "mostly driven by just 2-3 genes", run_go_bp_
## overlap.R). Applying the same logic to GSEA: filter to leadingEdge size
## >= MIN_LEADING_EDGE before calling a pathway significant for this figure
## (note this is separate from fgsea's minSize=15, which bounds the TESTED
## gene set size, not how many of those genes end up in the leading edge).
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_domain_split"
PADJ_THRESHOLD <- 0.25  # GSEA's own exploratory FDR q threshold (Subramanian et al. 2005), matching
                         # the project's other direct-GSEA/pairwise-GSEA figures -- more lenient than
                         # the padj<0.05 used for DE/ORA, deliberately, since this contrast (n=11 vs
                         # n=11) is underpowered relative to the vs-WT contrasts (n=11 vs n=539).
MIN_LEADING_EDGE <- 5
MAX_SHOW <- 40  # cap display for legibility (see 23_top10_paired_barplots.py's identical pattern);
                 # full retained set is always in fgsea_T1_vs_S1S4_combined_solid_sig.csv regardless

## =============================================================================
## A. Load the 4 per-collection direct T1_vs_S1S4 results (stat-ranked)
## =============================================================================
collections <- c("hallmark", "kegg", "reactome", "oncogenic")

all_results <- list()
for (coll in collections) {
  out_file <- file.path(RESULTS_DIR, "tables", ANALYSIS,
                         if (coll == "hallmark") "fgsea_T1_vs_S1S4_domain.csv" else paste0("fgsea_T1_vs_S1S4_", coll, ".csv"))
  if (!file.exists(out_file)) {
    stop("Missing ", out_file, " -- run run_gsea_master_statranked_fix.R first (it computes this contrast for all 4 collections).")
  }
  log_progress("Loading direct T1_vs_S1S4 ", coll, " result: ", out_file)
  res <- read.csv(out_file)
  res$leadingEdge <- strsplit(res$leadingEdge, ",\\s*")
  res$collection <- coll
  res$n_leading_edge <- lengths(res$leadingEdge)
  all_results[[coll]] <- res
}

pooled <- bind_rows(all_results)
log_progress("Pooled direct T1_vs_S1S4 GSEA: ", nrow(pooled), " pathway-tests across 4 collections")

## =============================================================================
## B. Significant AND not thin-evidence (leading edge >= MIN_LEADING_EDGE)
## =============================================================================
sig_thin  <- pooled %>% filter(padj < PADJ_THRESHOLD, n_leading_edge < MIN_LEADING_EDGE)
sig_solid <- pooled %>% filter(padj < PADJ_THRESHOLD, n_leading_edge >= MIN_LEADING_EDGE) %>% arrange(padj)

log_progress("Significant at padj<", PADJ_THRESHOLD, ": ", nrow(sig_thin) + nrow(sig_solid), " total -- ",
             nrow(sig_thin), " excluded as thin evidence (<", MIN_LEADING_EDGE, "-gene leading edge), ",
             nrow(sig_solid), " retained")
log_progress(paste(capture.output(print(table(sig_solid$collection))), collapse = "\n"))

save_table(sig_solid %>% select(pathway, collection, NES, padj, n_leading_edge, leadingEdge) %>%
             mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ",")),
           "fgsea_T1_vs_S1S4_combined_solid_sig", ANALYSIS)
save_table(sig_thin %>% select(pathway, collection, NES, padj, n_leading_edge),
           "fgsea_T1_vs_S1S4_combined_thin_excluded", ANALYSIS)

## =============================================================================
## C. Bar chart -- same NES/Set2/black-outline conventions as the other
## pathway figures in this project
## =============================================================================
SET2 <- c(Hallmark = "#66c2a5", KEGG = "#fc8d62", Reactome = "#8da0cb", Oncogenic = "#e78ac3")

clean_label <- function(pathway, coll) {
  prefix <- c(hallmark = "HALLMARK_", kegg = "KEGG_", reactome = "REACTOME_", oncogenic = "")
  label <- sub(paste0("^", prefix[coll]), "", pathway)
  label <- gsub("_", " ", label)
  if (nchar(label) > 50) label <- paste0(substr(label, 1, 47), "...")
  label
}

make_direct_barplot <- function(data, out_name, collections_label, max_show = MAX_SHOW) {
  if (nrow(data) == 0) {
    log_progress(out_name, ": no pathways survive both the significance and leading-edge filters -- skipping figure")
    return(invisible(NULL))
  }
  n_solid <- nrow(data)
  title_suffix <- ""
  plot_src <- data
  if (n_solid > max_show) {
    log_progress(out_name, ": capping display at top ", max_show, " of ", n_solid, " by |NES|")
    plot_src <- data %>% arrange(desc(abs(NES))) %>% head(max_show)
    title_suffix <- paste0("\n(top ", max_show, " of ", n_solid, " by effect size)")
  }

  plot_df <- plot_src %>%
    rowwise() %>%
    mutate(label = clean_label(pathway, collection),
           collection_disp = c(hallmark = "Hallmark", kegg = "KEGG", reactome = "Reactome", oncogenic = "Oncogenic")[collection]) %>%
    ungroup() %>%
    arrange(NES) %>%
    mutate(label = factor(label, levels = unique(label)))

  fig_height <- max(4, 0.42 * nrow(plot_df) + 2)
  n_collections <- length(unique(plot_df$collection_disp))

  p <- ggplot(plot_df, aes(x = NES, y = label, fill = collection_disp)) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_text(aes(label = paste0("FDR q=", ifelse(padj < 0.001, formatC(padj, format = "e", digits = 1), sprintf("%.3f", padj)),
                                  "  (", n_leading_edge, " genes)"),
                   hjust = ifelse(NES >= 0, -0.05, 1.05)),
              size = 3.3, colour = "grey30") +
    scale_fill_manual(values = SET2, name = "Collection", guide = if (n_collections > 1) "legend" else "none") +
    scale_x_continuous(expand = expansion(mult = 0.45)) +
    labs(x = "NES (direct T1 vs S1S4 GSEA)", y = NULL,
         title = paste0("Pathways directly distinguishing T1 from S1S4 (FDR q<", PADJ_THRESHOLD, ")", title_suffix),
         subtitle = paste0("Direct T1-vs-S1S4 contrast (not inferred from two separate vs-WT tests), ", collections_label, ",\n",
                            "leading edge >=", MIN_LEADING_EDGE, " genes. ", nrow(sig_thin %>% filter(collection %in% unique(data$collection))),
                            " nominally significant pathway(s) excluded as\nthin evidence (<", MIN_LEADING_EDGE, " driving genes).")) +
    theme_kcna5(base_size = 15) +
    theme(legend.position = "top", plot.title = element_text(size = 15))

  save_fig(p, out_name, ANALYSIS, width = 12.5, height = fig_height)
}

make_direct_barplot(sig_solid, "fgsea_T1_vs_S1S4_combined_barplot",
                     "pooled across Hallmark+KEGG+Reactome+Oncogenic", max_show = MAX_SHOW)
make_direct_barplot(sig_solid %>% filter(collection == "hallmark"), "fgsea_T1_vs_S1S4_hallmark_barplot",
                     "Hallmark only", max_show = MAX_SHOW)

log_progress("Direct T1 vs S1S4 combined GSEA figure(s) complete.")
