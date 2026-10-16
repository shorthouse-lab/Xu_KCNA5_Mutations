## ---------------------------------------------------------------------------
## Full pathway-universe UMAP landscape: every term from all 4 collections
## (Hallmark, KEGG, Reactome, oncogenic signatures -- ~2264 terms total) is
## embedded into 2D by gene-overlap (Jaccard) similarity, exactly analogous
## to run_pca_umap_mutant_positioning.R's patient-level UMAP but for
## PATHWAYS instead of patients. AllMutant_vs_WT enrichment results are then
## overlaid on top of this fixed background -- so the figure shows not just
## which pathways are significant, but where they sit relative to the
## entire space of testable biology, and whether significant hits cluster
## tightly (a coherent programme) or scatter (a diffuse one).
##
## Different from run_pathway_landscape_map_combined.R's emapplot, which
## only lays out and clusters the SIGNIFICANT subset -- this embeds
## everything once as a fixed reference map, which is the more direct
## analogue of the patient UMAP the user is comparing this to.
## ---------------------------------------------------------------------------
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(uwot)
})
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_domain_split"
COLLECTIONS <- c(hallmark = "Hallmark", kegg = "KEGG", reactome = "Reactome", oncogenic = "Oncogenic")

## =============================================================================
## A. Load all 4 collections, build one combined, collection-prefixed gene set list
## =============================================================================
all_genesets <- list()
for (coll in names(COLLECTIONS)) {
  pathways <- readRDS(file.path(RAW_DIR, paste0(coll, "_genesets.rds")))
  label <- COLLECTIONS[[coll]]
  names(pathways) <- sapply(names(pathways), function(nm) {
    clean_name <- gsub(paste0("^", toupper(coll), "_"), "", nm)
    clean_name <- gsub("_", " ", clean_name)
    paste0("[", label, "] ", clean_name)
  })
  all_genesets <- c(all_genesets, pathways)
}
log_progress("Total pathway universe: ", length(all_genesets), " terms across ", length(COLLECTIONS), " collections")

## =============================================================================
## B. Binary term x gene sparse matrix -> Jaccard similarity -> distance
## =============================================================================
## Dedup genes WITHIN each term first -- some source gene sets contain
## duplicate gene symbols (discovered the hard way: sparseMatrix() sums
## repeated (term, gene) index pairs rather than treating them as one
## membership flag, which silently corrupts the Jaccard calculation --
## union can go negative and Jaccard can exceed 1 or hit Inf, which is
## exactly what broke the first attempt at this).
all_genesets <- lapply(all_genesets, unique)
all_genes <- sort(unique(unlist(all_genesets)))
log_progress("Unique genes across all terms: ", length(all_genes))

term_names <- names(all_genesets)
i_idx <- unlist(lapply(seq_along(all_genesets), function(i) rep(i, length(all_genesets[[i]]))))
j_idx <- match(unlist(all_genesets), all_genes)
bin_mat <- sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(length(term_names), length(all_genes)),
                         dimnames = list(term_names, all_genes))
log_progress("Binary term x gene matrix: ", nrow(bin_mat), " x ", ncol(bin_mat))

log_progress("Computing pairwise Jaccard similarity (sparse crossprod)...")
intersection_mat <- tcrossprod(bin_mat)  # term x term, intersection counts
sizes <- Matrix::rowSums(bin_mat)
union_mat <- outer(sizes, sizes, "+") - as.matrix(intersection_mat)
jaccard_mat <- as.matrix(intersection_mat) / union_mat
jaccard_mat[is.na(jaccard_mat)] <- 0
dist_mat <- 1 - jaccard_mat
diag(dist_mat) <- 0

## =============================================================================
## C. UMAP embedding of the full pathway universe (precomputed distance)
## =============================================================================
log_progress("Running UMAP on the full ", nrow(dist_mat), "x", ncol(dist_mat), " term distance matrix...")
set.seed(42)
umap_coords <- umap(as.dist(dist_mat), n_neighbors = 15, min_dist = 0.1)
coords_df <- data.frame(term = term_names, UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2])
coords_df$collection <- sub("^\\[([A-Za-z]+)\\].*", "\\1", coords_df$term)
save_table(coords_df, "pathway_umap_universe_coords", ANALYSIS)
log_progress("UMAP embedding complete.")

## =============================================================================
## D. Overlay AllMutant_vs_WT enrichment results (reuses the combined GSEA
## already run in run_pathway_landscape_map_combined.R -- same term naming
## convention, so this is a direct join, no recomputation)
## =============================================================================
enrich_path <- file.path(RESULTS_DIR, "tables", ANALYSIS, "gsea_clusterprofiler_AllMutant_vs_WT_combined.csv")
stopifnot(file.exists(enrich_path))
enrich <- read.csv(enrich_path) %>% transmute(term = Description, NES, padj = p.adjust, tested = TRUE)

plot_df <- coords_df %>% left_join(enrich, by = "term") %>%
  mutate(tested = coalesce(tested, FALSE), sig = !is.na(padj) & padj < 0.05)

log_progress("Terms in universe: ", nrow(plot_df), " | tested (passed minGSSize/maxGSSize filters): ", sum(plot_df$tested),
             " | significant (padj<0.05): ", sum(plot_df$sig, na.rm = TRUE))

p <- ggplot() +
  geom_point(data = plot_df %>% filter(!sig), aes(x = UMAP1, y = UMAP2), colour = "grey85", size = 1) +
  geom_point(data = plot_df %>% filter(sig), aes(x = UMAP1, y = UMAP2, colour = NES), size = 3) +
  scale_colour_gradient2(low = "#2a78d6", mid = "#dddddd", high = "#d9410d", midpoint = 0, name = "NES") +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "Pathway universe landscape: AllMutant vs WT",
       subtitle = paste0("Grey = all ", sum(!plot_df$sig), " non-significant terms (Hallmark+KEGG+Reactome+Oncogenic, by gene-overlap similarity)\n",
                          "Coloured = ", sum(plot_df$sig, na.rm = TRUE), " significant terms (padj<0.05)")) +
  theme_kcna5(base_size = 15) +
  theme(panel.grid = element_blank())
save_fig(p, "pathway_umap_landscape_AllMutant_vs_WT", ANALYSIS, width = 11, height = 9)

## Faceted-by-collection version, so each database's enriched terms are
## separately visible against the SAME fixed background layout (same
## coordinates every panel -- only the highlighted/grey split changes)
p_facet <- ggplot() +
  geom_point(data = plot_df, aes(x = UMAP1, y = UMAP2), colour = "grey88", size = 0.8) +
  geom_point(data = plot_df %>% filter(sig), aes(x = UMAP1, y = UMAP2, colour = NES), size = 2.6) +
  scale_colour_gradient2(low = "#2a78d6", mid = "#dddddd", high = "#d9410d", midpoint = 0, name = "NES") +
  facet_wrap(~collection, ncol = 2) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "Pathway universe landscape by collection: AllMutant vs WT",
       subtitle = "Same fixed background layout in every panel; coloured points = significant in that collection (padj<0.05)") +
  theme_kcna5(base_size = 14) +
  theme(panel.grid = element_blank())
save_fig(p_facet, "pathway_umap_landscape_AllMutant_vs_WT_facet", ANALYSIS, width = 11, height = 10)

## =============================================================================
## E. Second overlay, using PER-COLLECTION corrected significance (the
## trusted, manuscript-headline numbers: 6+1+21+28=56 total) instead of the
## single combined correction across all 1578 terms used in sections A-D
## (which is stricter and only keeps 37) -- restores the pathways that are
## significant within their own collection but didn't survive the stricter
## combined test. Also coloured by SOURCE COLLECTION rather than NES, so the
## database-of-origin for each point is visible directly rather than only
## inferable from spatial clustering.
## =============================================================================
clean_pathway_name <- function(raw_name, coll, label) {
  clean <- gsub(paste0("^", toupper(coll), "_"), "", raw_name)
  clean <- gsub("_", " ", clean)
  paste0("[", label, "] ", clean)
}

trusted_files <- c(hallmark = "fgsea_AllMutant_vs_WT_hallmark.csv", kegg = "fgsea_AllMutant_vs_WT_kegg.csv",
                    reactome = "fgsea_AllMutant_vs_WT_reactome.csv", oncogenic = "fgsea_AllMutant_vs_WT_oncogenic.csv")

## Set2 (ColorBrewer/seaborn qualitative palette), 4 categories
SET2 <- c("Hallmark" = "#66c2a5", "KEGG" = "#fc8d62", "Reactome" = "#8da0cb", "Oncogenic" = "#e78ac3")

get_trusted_sig <- function(fdr_threshold) {
  bind_rows(lapply(names(trusted_files), function(coll) {
    df <- read.csv(file.path(RESULTS_DIR, "tables", "analysis3_transcriptomic", trusted_files[[coll]]))
    df %>% filter(padj < fdr_threshold) %>%
      transmute(term = clean_pathway_name(pathway, coll, COLLECTIONS[[coll]]), NES, padj,
                collection = COLLECTIONS[[coll]])
  }))
}

## Shared plot builder: grey background, Set2-filled + black-outlined points
## for the significant set, no axis ticks/text (not meaningful for a UMAP
## embedding -- only relative position/distance carries information). Point
## size scales with -log10(padj) -- matches the interactive UMAP dashboard,
## where more significant pathways render as visibly larger points.
build_landscape_plot <- function(plot_df_trusted, n_matched, title, subtitle) {
  sig_df <- plot_df_trusted %>% filter(sig) %>%
    mutate(neglog10padj = -log10(pmax(padj, 1e-16)))

  ggplot() +
    geom_point(data = plot_df_trusted %>% filter(!sig), aes(x = UMAP1, y = UMAP2), colour = "grey88", size = 1) +
    geom_point(data = sig_df, aes(x = UMAP1, y = UMAP2, fill = collection_sig, size = neglog10padj),
               shape = 21, colour = "black", stroke = 0.6) +
    scale_fill_manual(values = SET2, name = "Collection") +
    scale_size_continuous(range = c(1.8, 5.5), name = expression(-log[10](FDR~q))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title, subtitle = subtitle) +
    theme_kcna5(base_size = 15) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
}

## =============================================================================
## E. padj<0.05 within each collection's own correction (the manuscript's
## main-text numbers: 6+1+21+28=56)
## =============================================================================
trusted_sig <- get_trusted_sig(0.05)
log_progress("padj<0.05 per-collection significant terms (trusted, manuscript numbers): ", nrow(trusted_sig),
             " total\n", paste(capture.output(print(table(trusted_sig$collection))), collapse = "\n"))

plot_df_trusted <- coords_df %>% left_join(trusted_sig %>% select(term, NES, padj, collection_sig = collection), by = "term") %>%
  mutate(sig = !is.na(padj))
n_matched <- sum(plot_df_trusted$sig)
log_progress("Matched to UMAP background by term name: ", n_matched, " of ", nrow(trusted_sig),
             " (", nrow(trusted_sig) - n_matched, " failed to match -- check naming if this is large)")

p_collection <- build_landscape_plot(
  plot_df_trusted, n_matched,
  title = "Pathway universe landscape: AllMutant vs WT, coloured by source collection",
  subtitle = paste0("Grey = all ", sum(!plot_df_trusted$sig), " non-significant terms. Coloured = ", n_matched,
                     " significant terms (padj<0.05 within each\ncollection's own correction -- the main-text numbers, not a combined correction across collections)"))
save_fig(p_collection, "pathway_umap_landscape_AllMutant_vs_WT_by_collection", ANALYSIS, width = 12, height = 9)
save_table(trusted_sig, "pathway_umap_trusted_significant_terms", ANALYSIS)

## =============================================================================
## F. FDR q<0.25: the standard GSEA (Subramanian et al. 2005) exploratory
## threshold, more lenient than the padj<0.05 used everywhere else in this
## project. GSEA's own documentation recommends FDR q<0.25 specifically
## because its permutation-based FDR estimate is conservative relative to a
## typical differential-expression FDR -- this is NOT the same claim as
## "significant", and should be labelled as exploratory/suggestive if used,
## not substituted for the padj<0.05 result everywhere else.
## =============================================================================
trusted_sig_fdr25 <- get_trusted_sig(0.25)
log_progress("FDR q<0.25 per-collection terms (GSEA exploratory threshold): ", nrow(trusted_sig_fdr25),
             " total\n", paste(capture.output(print(table(trusted_sig_fdr25$collection))), collapse = "\n"))

plot_df_fdr25 <- coords_df %>% left_join(trusted_sig_fdr25 %>% select(term, NES, padj, collection_sig = collection), by = "term") %>%
  mutate(sig = !is.na(padj))
n_matched_fdr25 <- sum(plot_df_fdr25$sig)
log_progress("Matched to UMAP background: ", n_matched_fdr25, " of ", nrow(trusted_sig_fdr25))

p_fdr25 <- build_landscape_plot(
  plot_df_fdr25, n_matched_fdr25,
  title = "Pathway universe landscape: AllMutant vs WT, GSEA exploratory threshold (FDR q<0.25)",
  subtitle = paste0("Grey = non-significant. Coloured = ", n_matched_fdr25,
                     " terms at FDR q<0.25 (Subramanian et al. 2005's recommended GSEA\nexploratory cutoff -- more lenient than the padj<0.05 used elsewhere in this project; treat as suggestive, not confirmed"))
save_fig(p_fdr25, "pathway_umap_landscape_AllMutant_vs_WT_FDR25", ANALYSIS, width = 12, height = 9)
save_table(trusted_sig_fdr25, "pathway_umap_trusted_terms_FDR25", ANALYSIS)

log_progress("Pathway UMAP landscape complete.")
