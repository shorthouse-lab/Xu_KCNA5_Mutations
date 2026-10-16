## ---------------------------------------------------------------------------
## Does naturally low KCNA5 expression (no mutation) phenocopy the
## transcriptional effect of KCNA5 mutation? DESeq2 contrast: KCNA5-low vs
## KCNA5-high (bottom vs top quartile TPM), restricted to WT (non-mutant)
## patients only -- same restriction already used for the "expression_low"
## survival-analysis flag in define_groups.R (bottom quartile KCNA5 TPM
## among WT), extended here with the complementary top-quartile group and
## run as an actual DE contrast rather than just a survival covariate.
##
## Restricting to WT patients is deliberate: among mutants, KCNA5 TPM could
## itself be confounded by the mutation (e.g. NMD of a truncating allele,
## or simply correlated with mutation status for other reasons) -- WT-only
## isolates a pure "expression level" axis, orthogonal to mutation status,
## so the comparison below (concordance with AllMutant_vs_WT) is a genuine
## independent test rather than circular.
##
## Contrast direction chosen as Low vs High (not High vs Low) specifically
## so it's directly comparable to AllMutant_vs_WT: if KCNA5 mutations act
## through loss of channel function/expression, genes down in AllMutant
## (relative to WT) should also be down in Low (relative to High).
##
## Same design/covariates as every other DESeq2 contrast here: msi_score +
## TMB + group, Wald test, BH-adjusted.
## ---------------------------------------------------------------------------
library(DESeq2)
library(BiocParallel)
library(SummarizedExperiment)
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_transcriptomic"
N_WORKERS <- 8
register(MulticoreParam(workers = N_WORKERS))
log_progress("Registered MulticoreParam with ", N_WORKERS, " workers")

groups <- readRDS(file.path(PROCESSED_DIR, "patient_groups.rds"))
clin <- readRDS(file.path(RAW_DIR, "clinical_cbioportal.rds"))
se <- readRDS(file.path(RAW_DIR, "rnaseq_se.rds"))

tumor_se <- se[, substr(colnames(se), 14, 15) == "01"]
sample_patient <- substr(colnames(tumor_se), 1, 12)
libsize <- colSums(assay(tumor_se, "unstranded"))
keep_sample <- !duplicated(sample_patient[order(-libsize)])
tumor_se <- tumor_se[, order(-libsize)][, keep_sample]
colnames(tumor_se) <- substr(colnames(tumor_se), 1, 12)

msi_lookup <- clin %>% transmute(patient_id = patientId, msi_score = as.numeric(MSI_SENSOR_SCORE))

## =============================================================================
## Define expression quartile groups among WT patients only
## =============================================================================
wt_expr <- groups %>% filter(group == "WT", !is.na(kcna5_tpm))
q1 <- quantile(wt_expr$kcna5_tpm, 0.25)
q3 <- quantile(wt_expr$kcna5_tpm, 0.75)
log_progress("WT KCNA5 TPM quartiles: Q1=", round(q1, 3), ", Q3=", round(q3, 3),
             " (n=", nrow(wt_expr), " WT patients with expression data)")

coldata <- wt_expr %>%
  mutate(expr_group = case_when(
    kcna5_tpm <= q1 ~ "Low",
    kcna5_tpm >= q3 ~ "High",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(expr_group)) %>%
  mutate(expr_group = factor(expr_group, levels = c("High", "Low"))) %>%
  left_join(msi_lookup, by = "patient_id") %>%
  filter(patient_id %in% colnames(tumor_se), !is.na(msi_score), !is.na(TMB))

log_progress("DESeq2 cohort (WT only, top/bottom KCNA5 TPM quartile): ", nrow(coldata), " patients (",
             paste(names(table(coldata$expr_group)), table(coldata$expr_group), sep = "=", collapse = ", "), ")")

counts_mat <- round(assay(tumor_se, "unstranded")[, coldata$patient_id])
rownames(coldata) <- coldata$patient_id

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata,
                               design = ~ msi_score + TMB + expr_group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
log_progress("Genes retained after low-count filter: ", nrow(dds))

log_progress("Running DESeq() with parallel=TRUE (", N_WORKERS, " workers)...")
dds <- DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = MulticoreParam(N_WORKERS))
log_progress("DESeq() complete.")

gene_symbols <- rowData(tumor_se)$gene_name[match(rownames(dds), rownames(tumor_se))]

res <- results(dds, contrast = c("expr_group", "Low", "High"), alpha = 0.05,
                parallel = TRUE, BPPARAM = MulticoreParam(N_WORKERS))
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  mutate(gene_symbol = gene_symbols[match(ensembl_id, rownames(dds))], contrast = "KCNA5_Low_vs_High_WT") %>%
  relocate(gene_symbol, .after = ensembl_id)

n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
log_progress("KCNA5_Low_vs_High_WT: ", n_sig, " genes padj<0.05 (of ", sum(!is.na(res_df$padj)), " tested)")
save_table(res_df, "deseq2_KCNA5_Low_vs_High_WT", ANALYSIS)

## =============================================================================
## Hallmark fgsea on this contrast, for direct pathway-level comparison.
## GSEA computation itself lives in run_gsea_master_statranked_fix.R
## (tie-free ranking by DESeq2's Wald `stat` column, covering this contrast
## across all 4 collections) -- run that script after this one; it reads
## deseq2_KCNA5_Low_vs_High_WT.csv (written just above) and writes
## fgsea_KCNA5_Low_vs_High_WT_hallmark.csv (and kegg/reactome/oncogenic).
## No fgsea call lives in this script, so there is no ranking logic here to
## go stale.
## =============================================================================

## =============================================================================
## Concordance with AllMutant_vs_WT: does low expression phenocopy mutation?
## =============================================================================
allmut <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "deseq2_AllMutant_vs_WT.csv"))
allmut_hallmark <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "fgsea_AllMutant_vs_WT_hallmark.csv"))

## Gene-level: Spearman correlation of log2FC, restricted to genes tested in
## both. A handful of gene symbols map to >1 Ensembl ID in each table (35-57
## duplicates out of ~27k, mostly pseudogene/paralog pairs) -- dedup by
## largest |log2FoldChange| per symbol before joining, same convention as
## the build_ranks() dedup used for fgsea elsewhere, so the join is a clean
## 1:1 rather than many-to-many.
dedup_by_symbol <- function(df, lfc_col) {
  df %>% filter(!is.na(.data[[lfc_col]]), !is.na(gene_symbol), gene_symbol != "") %>%
    group_by(gene_symbol) %>% slice_max(abs(.data[[lfc_col]]), n = 1, with_ties = FALSE) %>% ungroup()
}
gene_merge <- inner_join(
  dedup_by_symbol(res_df, "log2FoldChange") %>% transmute(gene_symbol, lfc_expr = log2FoldChange, padj_expr = padj),
  dedup_by_symbol(allmut, "log2FoldChange") %>% transmute(gene_symbol, lfc_mut = log2FoldChange, padj_mut = padj),
  by = "gene_symbol"
)
lfc_cor <- cor.test(gene_merge$lfc_expr, gene_merge$lfc_mut, method = "spearman")
log_progress("Gene-level log2FC concordance (KCNA5_Low_vs_High vs AllMutant_vs_WT), n=", nrow(gene_merge),
             " genes tested in both: Spearman rho=", round(lfc_cor$estimate, 3), ", p=", signif(lfc_cor$p.value, 3))

## Pathway-level: NES concordance for the Hallmark pathways significant in AllMutant_vs_WT.
## fg_expr_df (fgsea_KCNA5_Low_vs_High_WT_hallmark.csv) is produced by
## run_gsea_master_statranked_fix.R -- run that script before this one.
expr_hallmark_path <- file.path(RESULTS_DIR, "tables", ANALYSIS, "fgsea_KCNA5_Low_vs_High_WT_hallmark.csv")
if (!file.exists(expr_hallmark_path)) {
  stop("Missing ", expr_hallmark_path, " -- run run_gsea_master_statranked_fix.R first.")
}
fg_expr_df <- read.csv(expr_hallmark_path)

sig_mut_pathways <- allmut_hallmark %>% filter(padj < 0.05) %>% pull(pathway)
pathway_merge <- inner_join(
  fg_expr_df %>% filter(pathway %in% sig_mut_pathways) %>% transmute(pathway, NES_expr = NES, padj_expr = padj),
  allmut_hallmark %>% filter(pathway %in% sig_mut_pathways) %>% transmute(pathway, NES_mut = NES, padj_mut = padj),
  by = "pathway"
) %>% arrange(NES_mut)
save_table(pathway_merge, "kcna5_expression_vs_mutation_hallmark_concordance", ANALYSIS)
log_progress("Hallmark NES concordance, for pathways significant in AllMutant_vs_WT:\n",
             paste(capture.output(print(pathway_merge, digits = 3)), collapse = "\n"))

concordance_summary <- data.frame(
  n_genes_tested_both = nrow(gene_merge),
  spearman_rho_logfc = as.numeric(lfc_cor$estimate),
  spearman_p = lfc_cor$p.value,
  n_mut_sig_pathways = length(sig_mut_pathways),
  n_concordant_direction = sum(sign(pathway_merge$NES_expr) == sign(pathway_merge$NES_mut))
)
save_table(concordance_summary, "kcna5_expression_vs_mutation_concordance_summary", ANALYSIS)
log_progress("Concordance summary:\n", paste(capture.output(print(concordance_summary, digits = 3)), collapse = "\n"))

saveRDS(list(dds = dds, res_df = res_df, coldata = coldata, gene_symbols = gene_symbols, fgsea_hallmark = fg_expr_df),
        file.path(PROCESSED_DIR, "analysis3_kcna5_expression_deseq2_results.rds"))
log_progress("KCNA5 expression-group DESeq2 contrast + concordance analysis complete.")
