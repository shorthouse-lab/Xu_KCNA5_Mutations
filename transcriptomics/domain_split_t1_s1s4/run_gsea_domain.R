## ---------------------------------------------------------------------------
## fgsea (Hallmark, 10,000 permutations) for the domain-split contrasts
## (T1_vs_WT, S1S4_vs_WT, T1_vs_S1S4), mirroring run_gsea.R.
##
## Also positions the single S5S6_mutant patient (TCGA-AY-4071, A470T --
## DESeq2 itself excludes this group entirely, see run_deseq2_contrasts_domain.R)
## descriptively against the T1/S1S4/WT distributions, using per-sample
## ssGSEA Hallmark scores rather than a group-level DESeq2/fgsea test that a
## single sample cannot support.
## ---------------------------------------------------------------------------
library(GSVA)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_domain_split"
hallmark <- readRDS(file.path(RAW_DIR, "hallmark_genesets.rds"))

## =============================================================================
## A. Load the 3 domain-split contrasts' Hallmark fgsea results
## GSEA computation itself lives in run_gsea_master_statranked_fix.R (tie-free
## ranking by DESeq2's Wald `stat` column -- run that script first). This
## script only consumes its output; there is no fgsea call here to go stale.
## =============================================================================
OUT_T1 <- file.path(RESULTS_DIR, "tables", ANALYSIS, "fgsea_T1_vs_WT_domain.csv")
OUT_S1S4 <- file.path(RESULTS_DIR, "tables", ANALYSIS, "fgsea_S1S4_vs_WT_domain.csv")
OUT_T1S4 <- file.path(RESULTS_DIR, "tables", ANALYSIS, "fgsea_T1_vs_S1S4_domain.csv")
if (!all(file.exists(OUT_T1, OUT_S1S4, OUT_T1S4))) {
  stop("Missing fgsea_{T1_vs_WT,S1S4_vs_WT,T1_vs_S1S4}_domain.csv -- run run_gsea_master_statranked_fix.R first.")
}
fgsea_T1_WT   <- read.csv(OUT_T1) %>% mutate(contrast = "T1_vs_WT")
fgsea_S1S4_WT <- read.csv(OUT_S1S4) %>% mutate(contrast = "S1S4_vs_WT")
fgsea_T1_S1S4 <- read.csv(OUT_T1S4) %>% mutate(contrast = "T1_vs_S1S4")

all_fgsea <- bind_rows(fgsea_T1_WT, fgsea_S1S4_WT, fgsea_T1_S1S4) %>%
  mutate(contrast = factor(contrast, levels = c("T1_vs_WT", "S1S4_vs_WT", "T1_vs_S1S4")),
         pathway = gsub("HALLMARK_", "", pathway))

pathway_order <- all_fgsea %>% group_by(pathway) %>% summarise(m = max(abs(NES))) %>%
  arrange(m) %>% pull(pathway)
all_fgsea$pathway <- factor(all_fgsea$pathway, levels = pathway_order)

p_dot <- ggplot(all_fgsea, aes(x = contrast, y = pathway, size = -log10(pmax(padj, 1e-10)), colour = NES)) +
  geom_point() +
  scale_colour_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  labs(x = NULL, y = NULL, size = "-log10(padj)", colour = "NES",
       title = "Hallmark GSEA across KCNA5 structural-domain contrasts") +
  theme_kcna5(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_dot, "fgsea_dotplot_domain_contrasts", ANALYSIS, width = 8, height = 11)

## =============================================================================
## B. Descriptive positioning of the single S5S6_mutant patient
## =============================================================================
log_progress("--- Descriptive: S5S6_mutant (n=1) vs T1/S1S4/WT ssGSEA score distributions ---")

groups_all <- readRDS(file.path(PROCESSED_DIR, "patient_groups_domain.rds"))
se <- readRDS(file.path(RAW_DIR, "rnaseq_se.rds"))
tumor_se <- se[, substr(colnames(se), 14, 15) == "01"]
sample_patient <- substr(colnames(tumor_se), 1, 12)
libsize <- colSums(assay(tumor_se, "unstranded"))
keep_sample <- !duplicated(sample_patient[order(-libsize)])
tumor_se <- tumor_se[, order(-libsize)][, keep_sample]
colnames(tumor_se) <- substr(colnames(tumor_se), 1, 12)

score_patients <- groups_all$patient_id[groups_all$domain_group %in% c("T1_mutant", "S1S4_mutant", "S5S6_mutant", "WT")]
score_patients <- intersect(score_patients, colnames(tumor_se))

tpm <- assay(tumor_se, "tpm_unstrand")[, score_patients]
gene_symbols_tpm <- rowData(tumor_se)$gene_name
rownames(tpm) <- gene_symbols_tpm
tpm <- tpm[!is.na(rownames(tpm)) & rownames(tpm) != "" & !duplicated(rownames(tpm)), ]
log2tpm <- log2(tpm + 1)

log_progress("Running ssGSEA (Hallmark) on ", ncol(log2tpm), " samples (incl. the single S5S6 patient)...")
gsvapar <- ssgseaParam(as.matrix(log2tpm), hallmark)
ssgsea_scores <- gsva(gsvapar)

score_df <- as.data.frame(t(ssgsea_scores)) %>%
  tibble::rownames_to_column("patient_id") %>%
  left_join(groups_all %>% select(patient_id, domain_group), by = "patient_id")

## Top pathways by significance in either T1_vs_WT or S1S4_vs_WT (from the formal test above)
top_pathways <- all_fgsea %>% filter(contrast %in% c("T1_vs_WT", "S1S4_vs_WT")) %>%
  group_by(pathway) %>% summarise(min_padj = min(padj)) %>% arrange(min_padj) %>% slice_head(n = 10) %>% pull(pathway)
top_pathways_full <- paste0("HALLMARK_", top_pathways)

plot_df <- score_df %>%
  filter(domain_group %in% c("T1_mutant", "S1S4_mutant", "WT")) %>%
  select(domain_group, all_of(top_pathways_full)) %>%
  tidyr::pivot_longer(-domain_group, names_to = "pathway", values_to = "score") %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway))

s5s6_df <- score_df %>%
  filter(domain_group == "S5S6_mutant") %>%
  select(all_of(top_pathways_full)) %>%
  tidyr::pivot_longer(everything(), names_to = "pathway", values_to = "score") %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway))

p_s5s6 <- ggplot(plot_df, aes(x = domain_group, y = score, fill = domain_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  geom_hline(data = s5s6_df, aes(yintercept = score), colour = "firebrick", linewidth = 0.8) +
  facet_wrap(~pathway, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "ssGSEA score", title = "Top Hallmark pathways: T1/S1S4/WT distributions vs the single S5S6 patient",
       subtitle = "Red line = TCGA-AY-4071 (the only S5S6_mutant patient, A470T) -- descriptive only, not a group test") +
  theme_kcna5(base_size = 12) + theme(legend.position = "none")
save_fig(p_s5s6, "ssgsea_s5s6_descriptive_position", ANALYSIS, width = 9, height = 11)

save_table(score_df, "ssgsea_hallmark_scores_domain_groups", ANALYSIS)
log_progress("GSEA + S5S6 descriptive positioning complete.")
