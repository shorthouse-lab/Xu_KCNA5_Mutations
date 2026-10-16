## ---------------------------------------------------------------------------
## Analysis 3, re-run with structural-domain-specific KCNA5 groups (T1,
## S1-S4 voltage-sensor, WT) instead of the original T1-vs-non-T1-vs-WT
## split -- following the structural findings in scripts/07_structural/
## (region-level dN/dS, S4 gating-charge enrichment, pore-domain
## hydrophobicity) with the transcriptomic question: do T1-mutant and
## S1-S4-mutant tumours actually look different from WT, and from each
## other, at the RNA-seq level?
##
## Same design and methodology as run_deseq2_contrasts.R: raw counts,
## MSI_SENSOR_SCORE + TMB as continuous covariates (character->numeric
## coercion is mandatory, see that script's header for the 6-hour-runtime
## bug this avoids), Wald test, BH-adjusted.
##
## S5S6_mutant (pore domain) has only 1 patient in TCGA-COADREAD (see
## define_domain_groups.R) -- DESeq2 needs within-group replication to
## estimate dispersion, so it is EXCLUDED from the model entirely here,
## not just from post-hoc testing. Its single sample is scored
## descriptively against the fitted model in run_gsea_domain.R instead.
## ---------------------------------------------------------------------------
library(DESeq2)
library(BiocParallel)
library(SummarizedExperiment)
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_domain_split"
N_WORKERS <- 8
register(MulticoreParam(workers = N_WORKERS))
log_progress("Registered MulticoreParam with ", N_WORKERS, " workers")

groups <- readRDS(file.path(PROCESSED_DIR, "patient_groups_domain.rds"))
clin <- readRDS(file.path(RAW_DIR, "clinical_cbioportal.rds"))
se <- readRDS(file.path(RAW_DIR, "rnaseq_se.rds"))

tumor_se <- se[, substr(colnames(se), 14, 15) == "01"]
sample_patient <- substr(colnames(tumor_se), 1, 12)
libsize <- colSums(assay(tumor_se, "unstranded"))
keep_sample <- !duplicated(sample_patient[order(-libsize)])
tumor_se <- tumor_se[, order(-libsize)][, keep_sample]
colnames(tumor_se) <- substr(colnames(tumor_se), 1, 12)

msi_lookup <- clin %>% transmute(patient_id = patientId, msi_score = as.numeric(MSI_SENSOR_SCORE))

coldata <- groups %>%
  filter(domain_group %in% c("T1_mutant", "S1S4_mutant", "WT")) %>%
  left_join(msi_lookup, by = "patient_id") %>%
  filter(patient_id %in% colnames(tumor_se), !is.na(msi_score), !is.na(TMB)) %>%
  mutate(group = factor(domain_group, levels = c("WT", "S1S4_mutant", "T1_mutant")))

log_progress("DESeq2 cohort: ", nrow(coldata), " patients (",
             paste(names(table(coldata$group)), table(coldata$group), sep = "=", collapse = ", "), ")")

counts_mat <- round(assay(tumor_se, "unstranded")[, coldata$patient_id])
rownames(coldata) <- coldata$patient_id

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata,
                               design = ~ msi_score + TMB + group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
log_progress("Genes retained after low-count filter: ", nrow(dds))

log_progress("Running DESeq() with parallel=TRUE (", N_WORKERS, " workers)...")
dds <- DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = MulticoreParam(N_WORKERS))
log_progress("DESeq() complete.")

gene_symbols <- rowData(tumor_se)$gene_name[match(rownames(dds), rownames(tumor_se))]

get_contrast <- function(grp1, grp2, label) {
  res <- results(dds, contrast = c("group", grp1, grp2), alpha = 0.05, parallel = TRUE, BPPARAM = MulticoreParam(N_WORKERS))
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("ensembl_id") %>%
    mutate(gene_symbol = gene_symbols[match(ensembl_id, rownames(dds))], contrast = label) %>%
    relocate(gene_symbol, .after = ensembl_id)
  n_sig <- sum(df$padj < 0.05, na.rm = TRUE)
  log_progress(label, ": ", n_sig, " genes padj<0.05 (of ", sum(!is.na(df$padj)), " tested)")
  df
}

res_T1_vs_WT     <- get_contrast("T1_mutant", "WT", "T1_vs_WT")
res_S1S4_vs_WT   <- get_contrast("S1S4_mutant", "WT", "S1S4_vs_WT")
res_T1_vs_S1S4   <- get_contrast("T1_mutant", "S1S4_mutant", "T1_vs_S1S4")

save_table(res_T1_vs_WT, "deseq2_T1_vs_WT_domain", ANALYSIS)
save_table(res_S1S4_vs_WT, "deseq2_S1S4_vs_WT_domain", ANALYSIS)
save_table(res_T1_vs_S1S4, "deseq2_T1_vs_S1S4_domain", ANALYSIS)

n_sig_direct <- sum(res_T1_vs_S1S4$padj < 0.05, na.rm = TRUE)
log_progress("T1 vs S1S4 direct contrast: ", n_sig_direct, " significant genes.")

saveRDS(list(dds = dds, res_T1_vs_WT = res_T1_vs_WT, res_S1S4_vs_WT = res_S1S4_vs_WT,
             res_T1_vs_S1S4 = res_T1_vs_S1S4, coldata = coldata, gene_symbols = gene_symbols),
        file.path(PROCESSED_DIR, "analysis3_domain_deseq2_results.rds"))
log_progress("Domain-split DESeq2 contrasts complete.")
