## ---------------------------------------------------------------------------
## Pooled contrast: ALL KCNA5 missense mutants (T1 + Non-T1, n=28) vs WT
## (n=539) -- the RNA-seq counterpart to the pooled "any alteration" survival
## analysis (run_analysis2_survival_any_alteration.R), but missense-mutation
## only (no amplification/non-missense), since those alteration types have
## no counts-based expression contrast to speak of in the same way and were
## separately shown to be an outlier-driven, biologically incoherent group
## when tried (Other_domain_mutant, see conversation -- not reused here).
##
## Same design/covariates as run_deseq2_contrasts.R for direct comparability:
## msi_score (continuous) + TMB (continuous) + group, Wald test, BH-adjusted.
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

coldata <- groups %>%
  filter(group %in% c("T1_mutant", "Non_T1_mutant", "WT")) %>%
  mutate(all_mutant = group != "WT",
         group2 = factor(ifelse(all_mutant, "AllMutant", "WT"), levels = c("WT", "AllMutant"))) %>%
  left_join(msi_lookup, by = "patient_id") %>%
  filter(patient_id %in% colnames(tumor_se), !is.na(msi_score), !is.na(TMB))

log_progress("DESeq2 cohort: ", nrow(coldata), " patients (",
             paste(names(table(coldata$group2)), table(coldata$group2), sep = "=", collapse = ", "), ")")

counts_mat <- round(assay(tumor_se, "unstranded")[, coldata$patient_id])
rownames(coldata) <- coldata$patient_id

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata,
                               design = ~ msi_score + TMB + group2)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
log_progress("Genes retained after low-count filter: ", nrow(dds))

log_progress("Running DESeq() with parallel=TRUE (", N_WORKERS, " workers)...")
dds <- DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = MulticoreParam(N_WORKERS))
log_progress("DESeq() complete.")

gene_symbols <- rowData(tumor_se)$gene_name[match(rownames(dds), rownames(tumor_se))]

res <- results(dds, contrast = c("group2", "AllMutant", "WT"), alpha = 0.05,
                parallel = TRUE, BPPARAM = MulticoreParam(N_WORKERS))
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  mutate(gene_symbol = gene_symbols[match(ensembl_id, rownames(dds))], contrast = "AllMutant_vs_WT") %>%
  relocate(gene_symbol, .after = ensembl_id)

n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
log_progress("AllMutant_vs_WT: ", n_sig, " genes padj<0.05 (of ", sum(!is.na(res_df$padj)), " tested)")

save_table(res_df, "deseq2_AllMutant_vs_WT", ANALYSIS)
saveRDS(list(dds = dds, res_AllMutant_vs_WT = res_df, coldata = coldata, gene_symbols = gene_symbols),
        file.path(PROCESSED_DIR, "analysis3_all_mutant_deseq2_results.rds"))
log_progress("AllMutant vs WT DESeq2 contrast complete.")
