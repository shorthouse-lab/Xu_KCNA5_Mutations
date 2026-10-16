## ---------------------------------------------------------------------------
## Patient grouping, applied once and reused by every downstream analysis.
##
## Definitions (verbatim from task spec):
##   T1-domain mutant : >=1 somatic MISSENSE mutation in KCNA5, residues 128-216
##   Non-T1 mutant    : missense mutation(s) outside 128-216, no T1 mutation
##   WT               : no somatic mutation AND no amplification in KCNA5
##   Expression-low   : bottom quartile KCNA5 TPM among WT patients (survival only)
##
## Note: patients whose ONLY KCNA5 alteration is non-missense (silent,
## frameshift, splice) satisfy none of the three definitions as literally
## written (they have a mutation, so aren't WT; it isn't missense, so aren't
## T1/non-T1) -- they are reported separately and excluded from the 3-group
## comparisons rather than silently folded into WT or non-T1.
## ---------------------------------------------------------------------------
library(maftools)
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

dir.create(PROCESSED_DIR, showWarnings = FALSE)

maf <- readRDS(file.path(RAW_DIR, "maf.rds"))
cnv_kcna5 <- readRDS(file.path(RAW_DIR, "cnv_kcna5.rds"))
cbio_kcna5 <- readRDS(file.path(RAW_DIR, "cbioportal_kcna5.rds"))

## --- 1. All patients in the cohort (from MAF clinical sample list) --------
all_samples <- as.character(maf@clinical.data$Tumor_Sample_Barcode)
all_patients <- unique(substr(all_samples, 1, 12))
log_progress("Cohort: ", length(all_patients), " unique patients (from MAF)")

## --- 2. KCNA5 mutation calls, parsed for position + type -------------------
kcna5_muts <- subsetMaf(maf = maf, genes = KCNA5_SYMBOL, mafObj = FALSE)
kcna5_muts <- as.data.frame(kcna5_muts)
kcna5_muts$patient_id <- substr(kcna5_muts$Tumor_Sample_Barcode, 1, 12)

## Parse residue position from HGVSp_Short, e.g. "p.D166N" -> 166,
## "p.T569Pfs*15" -> 569. Regex grabs the leading digits after the ref AA.
kcna5_muts$protein_pos <- as.integer(
  sub("^p\\.[A-Za-z*]+([0-9]+).*$", "\\1", kcna5_muts$HGVSp_Short)
)
kcna5_muts$is_missense <- kcna5_muts$Variant_Classification == "Missense_Mutation"
kcna5_muts$is_T1 <- kcna5_muts$is_missense &
  kcna5_muts$protein_pos >= T1_DOMAIN_START & kcna5_muts$protein_pos <= T1_DOMAIN_END

save_table(
  kcna5_muts %>% select(patient_id, Tumor_Sample_Barcode, Variant_Classification,
                         HGVSp_Short, protein_pos, is_missense, is_T1),
  "kcna5_mutation_calls", "02_patient_groups"
)

## --- 3. KCNA5 amplification calls (cBioPortal GISTIC: 2 = amplified) ------
cna <- cbio_kcna5$cna[[1]]
amp_patients <- unique(substr(gsub("-01$", "", cna$sampleId[cna$value == 2]), 1, 12))
log_progress("KCNA5-amplified patients (GISTIC=2): ", length(amp_patients))

## --- 4. Per-patient mutation summary ---------------------------------------
per_patient <- kcna5_muts %>%
  group_by(patient_id) %>%
  summarise(
    has_T1_missense     = any(is_T1),
    has_nonT1_missense  = any(is_missense & !is_T1),
    has_other_mutation  = any(!is_missense),
    .groups = "drop"
  )

groups <- data.frame(patient_id = all_patients, stringsAsFactors = FALSE) %>%
  left_join(per_patient, by = "patient_id") %>%
  mutate(
    has_T1_missense    = coalesce(has_T1_missense, FALSE),
    has_nonT1_missense = coalesce(has_nonT1_missense, FALSE),
    has_other_mutation = coalesce(has_other_mutation, FALSE),
    kcna5_amplified    = patient_id %in% amp_patients,
    has_any_mutation   = has_T1_missense | has_nonT1_missense | has_other_mutation,
    group = case_when(
      has_T1_missense                                          ~ "T1_mutant",
      has_nonT1_missense                                       ~ "Non_T1_mutant",
      !has_any_mutation & !kcna5_amplified                     ~ "WT",
      TRUE                                                     ~ "Excluded"  # non-missense-only mut, and/or amplified w/o WT-qualifying mutation status
    )
  )

log_progress("Group counts (pre expression-low):")
print(table(groups$group))

n_excluded <- sum(groups$group == "Excluded")
if (n_excluded > 0) {
  log_progress("NOTE: ", n_excluded, " patient(s) excluded from the 3 primary groups: ",
               "KCNA5 alteration present (non-missense mutation and/or amplification) ",
               "but doesn't meet the literal T1/Non-T1/WT definitions. See kcna5_mutation_calls table.")
}

## --- 5. TMB (mutations / Mb, all genes, not just KCNA5) --------------------
tmb_per_patient <- as.data.frame(maf@data) %>%
  mutate(patient_id = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  count(patient_id, name = "n_mutations") %>%
  mutate(TMB = n_mutations / TMB_CAPTURE_SIZE_MB)

groups <- groups %>% left_join(tmb_per_patient, by = "patient_id") %>%
  mutate(n_mutations = coalesce(n_mutations, 0L), TMB = coalesce(TMB, 0))

log_progress("TMB summary (mut/Mb, capture size = ", TMB_CAPTURE_SIZE_MB, " Mb): ",
             "median=", round(median(groups$TMB), 2),
             " IQR=[", round(quantile(groups$TMB, .25), 2), ", ", round(quantile(groups$TMB, .75), 2), "]")

## --- 6. Expression-low (bottom quartile KCNA5 TPM among WT, survival only) -
rnaseq_path <- file.path(RAW_DIR, "rnaseq_se.rds")
if (file.exists(rnaseq_path)) {
  library(SummarizedExperiment)
  se <- readRDS(rnaseq_path)
  ## rownames are versioned Ensembl IDs (e.g. "ENSG00000130037.5")
  kcna5_row <- grep(paste0("^", KCNA5_ENSEMBL, "\\."), rownames(se), value = TRUE)[1]
  kcna5_tpm <- assay(se, "tpm_unstrand")[kcna5_row, ]
  names(kcna5_tpm) <- substr(colnames(se), 1, 12)
  ## collapse to one value per patient (max across any multiple samples)
  kcna5_tpm_patient <- tapply(kcna5_tpm, names(kcna5_tpm), max, na.rm = TRUE)

  ## as.vector() drops the dim attribute tapply leaves behind, which
  ## readr::write_csv otherwise rejects as a "matrix column"
  groups$kcna5_tpm <- as.vector(kcna5_tpm_patient[groups$patient_id])

  wt_tpm <- groups$kcna5_tpm[groups$group == "WT"]
  q1_threshold <- quantile(wt_tpm, 0.25, na.rm = TRUE)
  groups$expression_low <- !is.na(groups$kcna5_tpm) &
    groups$group == "WT" & groups$kcna5_tpm <= q1_threshold

  log_progress("Expression-low (bottom quartile KCNA5 TPM among WT, threshold=",
               round(q1_threshold, 3), " TPM): ", sum(groups$expression_low), " patients")
} else {
  log_progress("rnaseq_se.rds not yet available -- skipping expression-low assignment for now. ",
               "Re-run this script once the RNA-seq download completes.")
  groups$kcna5_tpm <- NA_real_
  groups$expression_low <- NA
}

saveRDS(groups, file.path(PROCESSED_DIR, "patient_groups.rds"))
save_table(groups, "patient_groups", "02_patient_groups")
log_progress("Saved data_processed/patient_groups.rds (", nrow(groups), " patients)")
