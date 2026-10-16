## ---------------------------------------------------------------------------
## Extends patient_groups.rds with structural-domain-specific KCNA5 mutant
## groups (T1, S1-S4 voltage-sensor, S5-P-S6 pore), following up on the
## structural analyses in scripts/07_structural/ which found these three
## domains have distinct mutational signatures (region-level dN/dS, S4
## gating-charge enrichment, pore-radius/hydrophobicity context). This asks
## whether those structural differences translate into differences in the
## TCGA-COADREAD cohort's own genomic/clinical covariates (TMB, CMS, driver
## co-occurrence) -- re-running Analysis 1 with these finer groups instead of
## the original T1-vs-non-T1-vs-WT split.
##
## Domain boundaries match the structural analyses exactly:
##   T1    128-216   S1-S4 (VSD)  247-416   S5-P-S6 (pore)  433-613
## Positions outside all three (disordered regions, S0 231-240, linkers)
## fall to "Other_domain_mutant", reported but not used in the primary
## group comparisons (parallel to how the original script separately
## reports "Excluded" rather than silently folding edge cases into WT).
##
## IMPORTANT SAMPLE-SIZE CAVEAT, discovered by this script, not assumed in
## advance: TCGA-COADREAD (n=567 patients) has far fewer KCNA5 mutations
## than the pooled COSMIC database used for the structural work (28 unique
## missense-mutant patients here vs. ~69 unique mutated positions there).
## Splitting further by domain leaves T1 (n=11) and S1-S4 (n=11) reasonably
## matched, but S5-P-S6 has only ONE mutant patient (TCGA-AY-4071, A470T --
## the same residue independently flagged by the COSMIC structural work).
## Downstream scripts must NOT run formal group-comparison statistics on a
## group of size 1; report it descriptively only.
## ---------------------------------------------------------------------------
library(maftools)
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

groups <- readRDS(file.path(PROCESSED_DIR, "patient_groups.rds"))
maf <- readRDS(file.path(RAW_DIR, "maf.rds"))

kcna5_muts <- subsetMaf(maf = maf, genes = KCNA5_SYMBOL, mafObj = FALSE)
kcna5_muts <- as.data.frame(kcna5_muts)
kcna5_muts$patient_id <- substr(kcna5_muts$Tumor_Sample_Barcode, 1, 12)
kcna5_muts$protein_pos <- as.integer(sub("^p\\.[A-Za-z*]+([0-9]+).*$", "\\1", kcna5_muts$HGVSp_Short))
kcna5_muts$is_missense <- kcna5_muts$Variant_Classification == "Missense_Mutation"

S1S4_START <- 247; S1S4_END <- 416
S5S6_START <- 433; S5S6_END <- 613

mm <- kcna5_muts %>% filter(is_missense) %>%
  mutate(domain = case_when(
    protein_pos >= T1_DOMAIN_START & protein_pos <= T1_DOMAIN_END ~ "T1",
    protein_pos >= S1S4_START & protein_pos <= S1S4_END ~ "S1S4",
    protein_pos >= S5S6_START & protein_pos <= S5S6_END ~ "S5S6",
    TRUE ~ "Other"
  ))

## Sanity check carried over from interactive verification: no patient's
## KCNA5 missense mutation(s) span more than one of these domains, so
## precedence order doesn't matter here -- but check again on every run
## rather than assume it stays true if the MAF is ever refreshed.
multi_domain <- mm %>% distinct(patient_id, domain) %>% count(patient_id) %>% filter(n > 1)
if (nrow(multi_domain) > 0) {
  log_progress("WARNING: ", nrow(multi_domain), " patient(s) have missense mutations in >1 structural domain -- ",
               "precedence order (T1 > S1S4 > S5S6 > Other) applied below. Patients: ",
               paste(multi_domain$patient_id, collapse = ", "))
}

per_patient_domain <- mm %>%
  group_by(patient_id) %>%
  summarise(
    has_T1    = any(domain == "T1"),
    has_S1S4  = any(domain == "S1S4"),
    has_S5S6  = any(domain == "S5S6"),
    has_Other = any(domain == "Other"),
    .groups = "drop"
  )

groups_domain <- groups %>%
  left_join(per_patient_domain, by = "patient_id") %>%
  mutate(across(c(has_T1, has_S1S4, has_S5S6, has_Other), ~ coalesce(., FALSE))) %>%
  mutate(domain_group = case_when(
    has_T1    ~ "T1_mutant",
    has_S1S4  ~ "S1S4_mutant",
    has_S5S6  ~ "S5S6_mutant",
    has_Other ~ "Other_domain_mutant",
    group == "WT" ~ "WT",
    TRUE ~ "Excluded"   # non-missense-only mutation and/or amplified, same edge case as original grouping
  ))

log_progress("Domain group counts:")
print(table(groups_domain$domain_group))
log_progress("(For reference, TMB capture size used elsewhere: ", TMB_CAPTURE_SIZE_MB, " Mb)")
log_progress("CAVEAT: S5S6_mutant has only ", sum(groups_domain$domain_group == "S5S6_mutant"),
             " patient(s) -- report descriptively only, do not run group-comparison statistics on it.")

saveRDS(groups_domain, file.path(PROCESSED_DIR, "patient_groups_domain.rds"))
save_table(groups_domain, "patient_groups_domain", "02_patient_groups")
log_progress("Saved data_processed/patient_groups_domain.rds (", nrow(groups_domain), " patients)")
