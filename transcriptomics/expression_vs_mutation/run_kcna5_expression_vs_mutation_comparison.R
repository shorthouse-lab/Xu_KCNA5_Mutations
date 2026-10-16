## ---------------------------------------------------------------------------
## Figures for the "does KCNA5 mutation phenocopy low KCNA5 expression?"
## comparison (run_deseq2_kcna5_expression_groups.R): genome-wide log2FC
## concordance is null (rho=-0.02), but the 6 Hallmark pathways significant
## in AllMutant_vs_WT split 3-concordant (KRAS_SIGNALING_DN, INFLAMMATORY_
## RESPONSE, COMPLEMENT) / 3-discordant (E2F_TARGETS, MYC_TARGETS_V1,
## G2M_CHECKPOINT -- UP in low-expression WT, DOWN in mutants). Same figure
## conventions as the T1-vs-S1S4 domain comparison (run_domain_signature_
## comparison.R's scatter, run_pathway_comparison_dumbbell.R's dumbbell),
## reused here for a different pair of contrasts so the two comparisons
## read consistently across the manuscript.
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_transcriptomic"

expr <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "deseq2_KCNA5_Low_vs_High_WT.csv"))
mut <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "deseq2_AllMutant_vs_WT.csv"))

dedup_by_symbol <- function(df) {
  df %>% filter(!is.na(log2FoldChange), !is.na(gene_symbol), gene_symbol != "") %>%
    group_by(gene_symbol) %>% slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) %>% ungroup()
}

merged <- inner_join(
  dedup_by_symbol(expr) %>% transmute(gene_symbol, lfc_expr = log2FoldChange, padj_expr = padj),
  dedup_by_symbol(mut) %>% transmute(gene_symbol, lfc_mut = log2FoldChange, padj_mut = padj),
  by = "gene_symbol"
) %>%
  mutate(sig_status = case_when(
    padj_expr < 0.05 & padj_mut < 0.05 ~ "Significant in both",
    padj_expr < 0.05 ~ "KCNA5-expression-specific",
    padj_mut < 0.05 ~ "Mutation-specific",
    TRUE ~ "Neither"
  ))

ct <- cor.test(merged$lfc_expr, merged$lfc_mut, method = "spearman")
log_progress("Genome-wide log2FC concordance (n=", nrow(merged), "): Spearman rho=",
             round(ct$estimate, 3), ", p=", signif(ct$p.value, 3))
print(table(merged$sig_status))

# =============================================================================
# A. Genome-wide log2FC scatter
# =============================================================================
p_scatter <- ggplot(merged, aes(x = lfc_expr, y = lfc_mut, colour = sig_status)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = merged %>% filter(sig_status == "Neither"), alpha = 0.12, size = 0.7) +
  geom_point(data = merged %>% filter(sig_status != "Neither"), alpha = 0.6, size = 1.4) +
  scale_colour_manual(values = c("Significant in both" = "#1b9e77",
                                  "KCNA5-expression-specific" = "#7570b3",
                                  "Mutation-specific" = "#d95f02",
                                  "Neither" = "grey70"), name = NULL) +
  labs(x = "log2FC, KCNA5-low vs KCNA5-high (WT patients only)",
       y = "log2FC, AllMutant vs WT",
       title = "KCNA5 mutation does not transcriptionally phenocopy low KCNA5 expression",
       subtitle = sprintf("Genome-wide log2FC concordance: Spearman rho=%.3f (n=%d genes tested in both)",
                           ct$estimate, nrow(merged))) +
  theme_kcna5(base_size = 16) +
  theme(legend.position = "top", legend.box = "horizontal")
save_fig(p_scatter, "kcna5_expression_vs_mutation_concordance_scatter", ANALYSIS, width = 9, height = 7.5)

# =============================================================================
# B. Pathway-level dumbbell: the 6 Hallmark pathways significant in
#    AllMutant_vs_WT, NES in both contrasts
# =============================================================================
pathway_df <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "kcna5_expression_vs_mutation_hallmark_concordance.csv")) %>%
  mutate(pathway_label = gsub("HALLMARK_", "", pathway), pathway_label = gsub("_", " ", pathway_label),
         concordant = sign(NES_expr) == sign(NES_mut)) %>%
  arrange(NES_mut) %>%
  mutate(pathway_label = factor(pathway_label, levels = pathway_label))

long <- bind_rows(
  pathway_df %>% transmute(pathway_label, NES = NES_expr, padj = padj_expr, contrast = "KCNA5-low vs -high (WT)"),
  pathway_df %>% transmute(pathway_label, NES = NES_mut, padj = padj_mut, contrast = "AllMutant vs WT")
) %>%
  mutate(contrast = factor(contrast, levels = c("KCNA5-low vs -high (WT)", "AllMutant vs WT")))

p_dumbbell <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_segment(data = pathway_df, aes(x = NES_expr, xend = NES_mut, y = pathway_label, yend = pathway_label,
                                       colour = concordant), linewidth = 0.8) +
  geom_point(data = long, aes(x = NES, y = pathway_label, shape = contrast), size = 3.4, colour = "black") +
  scale_colour_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02"), guide = "none") +
  scale_shape_manual(values = c("KCNA5-low vs -high (WT)" = 1, "AllMutant vs WT" = 16), name = NULL) +
  labs(x = "Normalised enrichment score (NES)", y = NULL,
       title = "Mutation vs low expression: concordant or not?",
       subtitle = "Hallmark pathways significant in AllMutant vs WT.\nOrange connector = opposite direction between contrasts; green = same direction") +
  theme_kcna5(base_size = 16) +
  theme(legend.position = "top")
save_fig(p_dumbbell, "kcna5_expression_vs_mutation_hallmark_dumbbell", ANALYSIS, width = 10, height = 5.5)

log_progress("KCNA5 expression vs mutation comparison figures complete.")
