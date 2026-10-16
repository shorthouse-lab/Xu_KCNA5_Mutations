## ---------------------------------------------------------------------------
## How similar are the T1-mutant and S1S4-mutant transcriptional signatures,
## genome-wide, not just gene-by-gene or pathway-by-pathway? The per-contrast
## DESeq2/GSEA results (run_deseq2_contrasts_domain.R, run_gsea_domain.R)
## answer "how many genes/pathways are significant in each contrast" but not
## "do the two mutation classes push the transcriptome in the same direction
## overall, just with different strength/power, or genuinely different
## directions" -- that's a genome-wide concordance question, answered here
## two ways:
##   A. Correlate log2FoldChange(T1_vs_WT) against log2FoldChange(S1S4_vs_WT)
##      across every gene tested in both -- the standard way to compare two
##      differential expression signatures directly.
##   B. Same idea at the pathway level: correlate Hallmark NES(T1_vs_WT) vs
##      NES(S1S4_vs_WT) across all 50 gene sets, and tabulate which
##      significant sets are shared vs domain-specific.
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_domain_split"
res_list <- readRDS(file.path(PROCESSED_DIR, "analysis3_domain_deseq2_results.rds"))

## =============================================================================
## A. Genome-wide log2FC concordance, T1_vs_WT vs S1S4_vs_WT
## =============================================================================
t1 <- res_list$res_T1_vs_WT %>% filter(!is.na(padj)) %>%
  transmute(ensembl_id, gene_symbol, lfc_T1 = log2FoldChange, padj_T1 = padj)
s1s4 <- res_list$res_S1S4_vs_WT %>% filter(!is.na(padj)) %>%
  transmute(ensembl_id, gene_symbol, lfc_S1S4 = log2FoldChange, padj_S1S4 = padj)

merged <- inner_join(t1, s1s4, by = "ensembl_id") %>%
  mutate(gene_symbol = coalesce(gene_symbol.x, gene_symbol.y)) %>%
  select(-gene_symbol.x, -gene_symbol.y) %>%
  mutate(sig_status = case_when(
    padj_T1 < 0.05 & padj_S1S4 < 0.05 ~ "Significant in both",
    padj_T1 < 0.05 ~ "T1-specific",
    padj_S1S4 < 0.05 ~ "S1S4-specific",
    TRUE ~ "Neither"
  ))

log_progress("Genes tested in both T1_vs_WT and S1S4_vs_WT: ", nrow(merged))

ct <- cor.test(merged$lfc_T1, merged$lfc_S1S4, method = "spearman")
log_progress("Genome-wide log2FC concordance (all ", nrow(merged), " genes tested in both): ",
             "Spearman rho=", round(ct$estimate, 3), ", p=", signif(ct$p.value, 3))

## restricted to genes significant in at least one contrast (where the
## direction actually matters biologically, not just noise around 0)
sig_either <- merged %>% filter(sig_status != "Neither")
ct_sig <- cor.test(sig_either$lfc_T1, sig_either$lfc_S1S4, method = "spearman")
log_progress("Concordance restricted to genes significant in EITHER contrast (n=", nrow(sig_either), "): ",
             "Spearman rho=", round(ct_sig$estimate, 3), ", p=", signif(ct_sig$p.value, 3))

save_table(merged, "domain_signature_gene_concordance", ANALYSIS)

overlap_counts <- table(merged$sig_status)
log_progress("Gene-level overlap:\n", paste(capture.output(print(overlap_counts)), collapse = "\n"))
save_table(as.data.frame(overlap_counts) %>% setNames(c("category", "n_genes")), "domain_signature_overlap_counts", ANALYSIS)

label_genes <- merged %>% filter(sig_status == "Significant in both") %>%
  bind_rows(merged %>% filter(sig_status != "Neither") %>% slice_max(abs(lfc_T1) + abs(lfc_S1S4), n = 8))

p_scatter <- ggplot(merged, aes(x = lfc_T1, y = lfc_S1S4, colour = sig_status)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = merged %>% filter(sig_status == "Neither"), alpha = 0.15, size = 0.8) +
  geom_point(data = merged %>% filter(sig_status != "Neither"), alpha = 0.7, size = 1.6) +
  geom_text_repel(data = label_genes, aes(label = gene_symbol), size = 2.8, max.overlaps = 20) +
  scale_colour_manual(values = c("Neither" = "grey70", "T1-specific" = "steelblue",
                                  "S1S4-specific" = "darkorange", "Significant in both" = "firebrick")) +
  labs(x = "log2FC, T1-mutant vs WT", y = "log2FC, S1S4-mutant vs WT", colour = NULL,
       title = "Genome-wide transcriptional concordance: T1 vs S1S4 (both relative to WT)",
       subtitle = paste0("Spearman rho=", round(ct$estimate, 3), " (all genes tested in both, n=", nrow(merged), "), p=", signif(ct$p.value, 3))) +
  theme_kcna5(base_size = 16)
save_fig(p_scatter, "domain_signature_concordance_scatter", ANALYSIS, width = 9, height = 7)

## =============================================================================
## B. Pathway-level (Hallmark NES) concordance
## =============================================================================
fgsea_T1   <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "fgsea_T1_vs_WT_domain.csv"))
fgsea_S1S4 <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "fgsea_S1S4_vs_WT_domain.csv"))

## FDR q<0.25 -- GSEA's own exploratory threshold (Subramanian et al. 2005),
## used throughout this project's other direct/pairwise GSEA figures rather
## than the padj<0.05 used for DE/ORA (see run_pathway_umap_landscape.R,
## run_gsea_t1_vs_s1s4_direct_combined.R). More lenient than 0.05, matches
## the convention this figure is being brought in line with.
FDR_THRESHOLD <- 0.25

nes_merged <- inner_join(
  fgsea_T1 %>% transmute(pathway, NES_T1 = NES, padj_T1 = padj),
  fgsea_S1S4 %>% transmute(pathway, NES_S1S4 = NES, padj_S1S4 = padj),
  by = "pathway"
) %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         sig_status = case_when(
           padj_T1 < FDR_THRESHOLD & padj_S1S4 < FDR_THRESHOLD ~ "Significant in both",
           padj_T1 < FDR_THRESHOLD ~ "T1-specific",
           padj_S1S4 < FDR_THRESHOLD ~ "S1S4-specific",
           TRUE ~ "Neither"
         ))

ct_nes <- cor.test(nes_merged$NES_T1, nes_merged$NES_S1S4, method = "spearman")
log_progress("Pathway-level (Hallmark NES) concordance across all ", nrow(nes_merged), " sets: ",
             "Spearman rho=", round(ct_nes$estimate, 3), ", p=", signif(ct_nes$p.value, 3))
save_table(nes_merged, "domain_signature_pathway_concordance", ANALYSIS)

## Same T1/S1S4 palette as the venn diagrams (venn_genes/pathways_*_T1_vs_S1S4,
## Dark2 orange/purple) for consistency across every T1-vs-S1S4 figure.
## "Significant in both" uses Dark2 green (#1b9e77) -- the third colour in
## the same ColorBrewer Dark2 set the venn's orange/purple come from, so it
## reads as "this palette family" rather than an unrelated colour.
DOMAIN_COLOURS <- c("Neither" = "grey70", "T1-specific" = "#d95f02",
                     "S1S4-specific" = "#7570b3", "Significant in both" = "#1b9e77")

## Labels only on the most informative points: both pathways shared between
## contrasts, every S1S4-specific pathway (there are few), and just the
## top/bottom 10 by NES for the much larger T1-specific set -- avoids the
## label clutter of annotating all 19 T1-specific points.
t1_specific <- nes_merged %>% filter(sig_status == "T1-specific")
label_df <- bind_rows(
  nes_merged %>% filter(sig_status == "Significant in both"),
  nes_merged %>% filter(sig_status == "S1S4-specific"),
  t1_specific %>% slice_max(NES_T1, n = 10),
  t1_specific %>% slice_min(NES_T1, n = 10)
) %>% distinct(pathway, .keep_all = TRUE)

nes_lim <- max(abs(c(nes_merged$NES_T1, nes_merged$NES_S1S4))) * 1.1

p_nes <- ggplot(nes_merged, aes(x = NES_T1, y = NES_S1S4, colour = sig_status)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2) +
  geom_text_repel(data = label_df, aes(label = pathway), size = 2.8, max.overlaps = 30) +
  scale_colour_manual(values = DOMAIN_COLOURS, name = NULL) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(x = "NES, T1-mutant vs WT", y = "NES, S1S4-mutant vs WT", colour = NULL,
       title = "Hallmark pathway-level concordance: T1 vs S1S4 (both relative to WT)",
       subtitle = paste0("Spearman rho=", round(ct_nes$estimate, 3), " (all 50 sets), p=", signif(ct_nes$p.value, 3),
                          "  |  coloured at FDR q<", FDR_THRESHOLD)) +
  theme_kcna5(base_size = 16)
save_fig(p_nes, "domain_signature_concordance_pathways", ANALYSIS, width = 9, height = 9)

pathway_overlap <- table(nes_merged$sig_status)
log_progress("Pathway-level overlap:\n", paste(capture.output(print(pathway_overlap)), collapse = "\n"))
log_progress("Shared significant pathways: ", paste(nes_merged$pathway[nes_merged$sig_status == "Significant in both"], collapse = ", "))

## =============================================================================
## C. Same concordance plot, pooled across all 4 GSEA collections (Hallmark +
## KEGG + Reactome + Oncogenic) instead of Hallmark alone -- same FDR q<0.25
## threshold, same domain palette/labelling scheme, reused as a function so
## the two figures stay visually identical apart from scope.
## =============================================================================
make_concordance_plot <- function(data, title, out_name, label_cap = Inf) {
  ## Labels: top/bottom `label_cap` by NES_T1 within EACH category (not just
  ## T1-specific) -- with the full pooled universe, "Significant in both"
  ## (134 pts) and "S1S4-specific" (70 pts) are themselves too large to
  ## label in full, so the same most-extreme-points-only logic used for
  ## T1-specific in the Hallmark-only figure is applied to all three
  ## coloured categories here. label_cap = Inf (the Hallmark-only default)
  ## keeps the original "label everything in both/S1S4-specific" behaviour
  ## for that smaller figure.
  cap_group <- function(df) {
    if (nrow(df) <= 2 * label_cap) return(df)
    bind_rows(df %>% slice_max(NES_T1, n = label_cap), df %>% slice_min(NES_T1, n = label_cap))
  }
  t1_spec <- data %>% filter(sig_status == "T1-specific")
  label_df <- bind_rows(
    cap_group(data %>% filter(sig_status == "Significant in both")),
    cap_group(data %>% filter(sig_status == "S1S4-specific")),
    t1_spec %>% slice_max(NES_T1, n = 10),
    t1_spec %>% slice_min(NES_T1, n = 10)
  ) %>% distinct(pathway_id, .keep_all = TRUE)

  lim <- max(abs(c(data$NES_T1, data$NES_S1S4))) * 1.1
  ct <- cor.test(data$NES_T1, data$NES_S1S4, method = "spearman")

  p <- ggplot(data, aes(x = NES_T1, y = NES_S1S4, colour = sig_status)) +
    geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = 2) +
    geom_text_repel(data = label_df, aes(label = label), size = 2.6, max.overlaps = 40) +
    scale_colour_manual(values = DOMAIN_COLOURS, name = NULL) +
    coord_fixed(ratio = 1, xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = "NES, T1-mutant vs WT", y = "NES, S1S4-mutant vs WT", colour = NULL,
         title = title,
         subtitle = paste0("Spearman rho=", round(ct$estimate, 3), " (all ", nrow(data), " sets), p=", signif(ct$p.value, 3),
                            "  |  coloured at FDR q<", FDR_THRESHOLD)) +
    theme_kcna5(base_size = 15) +
    theme(plot.title = element_text(size = 14))
  save_fig(p, out_name, ANALYSIS, width = 10.5, height = 9)
  list(data = data, cor = ct)
}

COLLECTION_FILES <- list(
  hallmark  = list(prefix = "HALLMARK_",  t1 = "fgsea_T1_vs_WT_domain.csv",   s1s4 = "fgsea_S1S4_vs_WT_domain.csv"),
  kegg      = list(prefix = "KEGG_",      t1 = "fgsea_T1_vs_WT_kegg.csv",     s1s4 = "fgsea_S1S4_vs_WT_kegg.csv"),
  reactome  = list(prefix = "REACTOME_",  t1 = "fgsea_T1_vs_WT_reactome.csv", s1s4 = "fgsea_S1S4_vs_WT_reactome.csv"),
  oncogenic = list(prefix = "",           t1 = "fgsea_T1_vs_WT_oncogenic.csv",s1s4 = "fgsea_S1S4_vs_WT_oncogenic.csv")
)

pooled <- bind_rows(lapply(names(COLLECTION_FILES), function(coll) {
  spec <- COLLECTION_FILES[[coll]]
  t1_c   <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, spec$t1))
  s1s4_c <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, spec$s1s4))
  inner_join(
    t1_c %>% transmute(pathway, NES_T1 = NES, padj_T1 = padj),
    s1s4_c %>% transmute(pathway, NES_S1S4 = NES, padj_S1S4 = padj),
    by = "pathway"
  ) %>%
    mutate(collection = coll, pathway_id = paste0(coll, ":", pathway),
           label = if (nchar(spec$prefix) > 0) sub(paste0("^", spec$prefix), "", pathway) else pathway,
           label = gsub("_", " ", label))
}))

pooled <- pooled %>%
  mutate(sig_status = case_when(
    padj_T1 < FDR_THRESHOLD & padj_S1S4 < FDR_THRESHOLD ~ "Significant in both",
    padj_T1 < FDR_THRESHOLD ~ "T1-specific",
    padj_S1S4 < FDR_THRESHOLD ~ "S1S4-specific",
    TRUE ~ "Neither"
  ))

log_progress("Pooled 4-collection pathway universe: ", nrow(pooled), " pathway-tests")
ct_pooled <- make_concordance_plot(
  pooled, "Pathway-level concordance: T1 vs S1S4 (4 collections pooled)",
  "domain_signature_concordance_pathways_all_collections", label_cap = 5
)$cor
save_table(pooled, "domain_signature_pathway_concordance_all_collections", ANALYSIS)

pooled_overlap <- table(pooled$sig_status)
log_progress("Pooled pathway-level overlap:\n", paste(capture.output(print(pooled_overlap)), collapse = "\n"))
log_progress("Pooled Spearman rho=", round(ct_pooled$estimate, 3), ", p=", signif(ct_pooled$p.value, 3))

saveRDS(list(gene_concordance = merged, gene_cor = ct, gene_cor_sig_only = ct_sig,
             pathway_concordance = nes_merged, pathway_cor = ct_nes,
             pathway_concordance_pooled = pooled, pathway_cor_pooled = ct_pooled),
        file.path(PROCESSED_DIR, "analysis3_domain_signature_comparison_results.rds"))
log_progress("Domain signature comparison complete.")
