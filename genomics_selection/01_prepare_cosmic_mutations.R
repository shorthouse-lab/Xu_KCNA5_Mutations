## ---------------------------------------------------------------------------
## Regenerate the whole-protein KCNA5 missense mutation list from raw COSMIC
## data (not checked into the Xu_KCNA5_Mutations-main repo due to licensing).
## Replicates the exact filtering/dedup logic from
## Xu_KCNA5_Mutations-main/Clustering/nmc.R:
##   - MUTATION_DESCRIPTION == "missense_variant" & GENE_SYMBOL == "KCNA5"
##   - distinct(SAMPLE_NAME, GENOMIC_MUTATION_ID, .keep_all = TRUE)
##   - position parsed from MUTATION_AA (e.g. "p.T417M" -> 417)
##
## Source file: /Users/davidshorthouse/Documents/Data/COSMIC/
##   Large_intestine_COSMIC_mutations.csv -- already tissue-filtered to large
##   intestine (matches the Xu repo's "LI" naming); histology
##   (carcinoma/adenocarcinoma) filtering could not be independently verified
##   since this export has no histology column, unlike what the Xu repo's
##   README describes for their own COSMIC download. Documented as an
##   assumption, not a verified match.
## ---------------------------------------------------------------------------
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
RAW_COSMIC <- "/tmp/kcna5_cosmic_raw.csv"  # pre-filtered via awk (GENE_SYMBOL==KCNA5) from the 1.2GB source

raw <- read.csv(RAW_COSMIC, stringsAsFactors = FALSE)
log_progress("Raw KCNA5 COSMIC records (large intestine): ", nrow(raw))
log_progress("Columns: ", paste(colnames(raw), collapse = ", "))

gene_data <- raw %>%
  filter(MUTATION_DESCRIPTION == "missense_variant", GENE_SYMBOL == "KCNA5") %>%
  distinct(SAMPLE_NAME, GENOMIC_MUTATION_ID, .keep_all = TRUE) %>%
  select(SAMPLE_NAME, GENOMIC_MUTATION_ID, MUTATION_AA)

## Parse position from "p.T417M" -> strip "p." (substring from char 3), then
## strip first (WT aa) and last (mut aa) characters -- identical to nmc.R
gene_data$aa_change <- substring(gene_data$MUTATION_AA, 3)
gene_data$position <- as.integer(substr(gene_data$aa_change, 2, nchar(gene_data$aa_change) - 1))
gene_data$wt_aa <- substr(gene_data$aa_change, 1, 1)
gene_data$mut_aa <- substr(gene_data$aa_change, nchar(gene_data$aa_change), nchar(gene_data$aa_change))

log_progress("Whole-protein KCNA5 missense mutations after dedup: ", nrow(gene_data),
             " (Xu repo / task prompt reference figure: 103 -- see script header for why an exact ",
             "match isn't guaranteed; reporting the actual regenerated count transparently)")

save_table(gene_data, "cosmic_kcna5_whole_protein_mutations", ANALYSIS)

## --- T1-domain subset (128-216, per explicit task instruction -- NOT the ---
## --- Xu repo's own internal 121-216 convention used in their T1 file)    ---
t1_data <- gene_data %>% filter(position >= T1_DOMAIN_START, position <= T1_DOMAIN_END)
log_progress("T1-domain (", T1_DOMAIN_START, "-", T1_DOMAIN_END, ") missense mutations: ", nrow(t1_data),
             " at ", length(unique(t1_data$position)), " unique positions: ",
             paste(sort(unique(t1_data$position)), collapse = ", "))
save_table(t1_data, "cosmic_kcna5_T1_domain_mutations", ANALYSIS)

saveRDS(list(whole_protein = gene_data, t1_domain = t1_data),
        file.path(PROCESSED_DIR, "structural_cosmic_mutations.rds"))
log_progress("COSMIC mutation list preparation complete.")
