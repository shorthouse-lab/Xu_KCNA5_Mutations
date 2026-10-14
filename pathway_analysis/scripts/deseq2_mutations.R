# =============================================================================
# deseq2_mutations.R
# BIOL0041 - KCNA5 Pathway Analysis in TCGA-COAD
#
# Purpose: DESeq2 differential expression analysis
#          1. All KCNA5 mutant vs wild-type tumours
#          2. T1 domain mutant (residues 128-216) vs wild-type tumours
# Input:   data/cleaned/TCGA-COAD_counts_cleaned.rds
#          data/cleaned/tumour_samples.rds
#          data/raw/KCNA5_mutations.txt
# Output:  results/deseq2/DESeq2_KCNA5_allMutant_vs_WT.csv
#          results/deseq2/DESeq2_KCNA5_T1mutant_vs_WT.csv
#          results/deseq2/dds_allMutant_vs_WT.rds
#          results/deseq2/dds_T1mutant_vs_WT.rds
# =============================================================================

library(DESeq2)
library(dplyr)
library(tibble)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
counts_raw     <- readRDS("data/cleaned/TCGA-COAD_counts_cleaned.rds")
tumour_samples <- readRDS("data/cleaned/tumour_samples.rds")

cat("Counts matrix loaded:", dim(counts_raw), "\n")
cat("Tumour samples:", length(tumour_samples), "\n")

# -----------------------------------------------------------------------------
# Load and clean KCNA5 mutation status
# -----------------------------------------------------------------------------
kcna5_status <- read.delim("data/raw/KCNA5_mutations.txt",
                             check.names = FALSE,
                             header      = TRUE)

kcna5_status_clean <- kcna5_status %>%
    filter(KCNA5 != "NP")

cat("Mutation status after removing NP:\n")
print(table(kcna5_status_clean$KCNA5))

# -----------------------------------------------------------------------------
# Identify mutant and WT samples
# -----------------------------------------------------------------------------
mutant_all <- paste0(
    kcna5_status_clean$SAMPLE_ID[kcna5_status_clean$KCNA5 != "WT"], "A")
wt_all     <- paste0(
    kcna5_status_clean$SAMPLE_ID[kcna5_status_clean$KCNA5 == "WT"], "A")

mutant_in_counts <- mutant_all[mutant_all %in% colnames(counts_raw) &
                                mutant_all %in% tumour_samples]
wt_in_counts     <- wt_all[wt_all %in% colnames(counts_raw) &
                             wt_all %in% tumour_samples]

cat("\nAll mutant tumour samples:", length(mutant_in_counts), "\n")
cat("WT tumour samples:", length(wt_in_counts), "\n")

# Save WT samples for use in other scripts
saveRDS(wt_in_counts, "data/cleaned/wt_tumour_samples.rds")

# -----------------------------------------------------------------------------
# Identify T1 domain mutations (residues 120-216)
# Restricted to structurally resolved region
# Residues 1-127 absent from crystal structure — excluded
# -----------------------------------------------------------------------------
t1_mutations <- data.frame(
    sample   = mutant_in_counts,
    mutation = kcna5_status_clean$KCNA5[
        match(gsub("A$", "", mutant_in_counts),
              kcna5_status_clean$SAMPLE_ID)]
) %>%
    mutate(
        residue = as.numeric(regmatches(mutation,
                                         regexpr("[0-9]+", mutation)))
    )

t1_samples        <- t1_mutations %>%
    filter(residue >= 120 & residue <= 216)
t1_mutant_samples <- t1_samples$sample

cat("\nT1 domain mutant samples (residues 120-216):\n")
print(t1_samples)
cat("Number of T1 domain mutant samples:", nrow(t1_samples), "\n")

# Binomial test - are T1 domain mutations enriched?
binom_result <- binom.test(
    x           = nrow(t1_samples),
    n           = length(mutant_in_counts),
    p           = 89/613,
    alternative = "greater"
)
cat("\nBinomial test for T1 domain mutation enrichment:\n")
cat("Observed:", nrow(t1_samples), "/", length(mutant_in_counts),
    "(", round(100 * nrow(t1_samples)/length(mutant_in_counts), 1), "%)\n")
cat("Expected under random distribution: 14.5%\n")
cat("p-value:", binom_result$p.value, "\n")

# Save sample lists
saveRDS(mutant_in_counts, "data/cleaned/all_mutant_samples.rds")
saveRDS(t1_mutant_samples, "data/cleaned/t1_mutant_samples.rds")
saveRDS(t1_samples, "data/cleaned/t1_samples_df.rds")

# -----------------------------------------------------------------------------
# Helper function: run DESeq2 and extract results
# -----------------------------------------------------------------------------
run_deseq2 <- function(counts_raw, 
                        case_samples, 
                        control_samples,
                        case_label    = "Mutant",
                        control_label = "WT",
                        output_prefix = "results/deseq2/DESeq2") {
    
    cat("\n=== DESeq2:", case_label, "vs", control_label, "===\n")
    
    # Subset counts
    counts_subset <- counts_raw[, c(case_samples, control_samples)]
    
    # Build metadata
    metadata <- data.frame(
        sample    = c(case_samples, control_samples),
        condition = factor(c(rep(case_label,    length(case_samples)),
                             rep(control_label, length(control_samples))),
                           levels = c(control_label, case_label))
    )
    rownames(metadata) <- metadata$sample
    
    cat("Sample counts:\n")
    print(table(metadata$condition))
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
        countData = counts_subset,
        colData   = metadata,
        design    = ~condition
    )
    
    # Filter lowly expressed genes
    keep <- rowSums(counts(dds) >= 10) >= 10
    dds  <- dds[keep, ]
    cat("Genes retained after filtering:", nrow(dds), "\n")
    
    # Run DESeq2
    cat("Running DESeq2...\n")
    dds <- DESeq(dds)
    
    # Extract results
    res <- results(dds,
                   contrast = c("condition", case_label, control_label),
                   alpha    = 0.05)
    
    cat("\nDESeq2 results summary:\n")
    summary(res)
    
    # Convert to dataframe
    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene") %>%
        filter(!is.na(padj)) %>%
        arrange(padj)
    
    cat("\nSignificant genes (padj < 0.05):",
        sum(res_df$padj < 0.05), "\n")
    cat("Upregulated in", case_label, ":",
        sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0), "\n")
    cat("Downregulated in", case_label, ":",
        sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0), "\n")
    
    # Check genes of interest
    genes_of_interest <- c("KCNA5", "FZD1", "FZD2", "FZD6", "FZD8",
                            "CTNNB1", "AXIN1", "MYC", "NOTCH4", "GNAI1")
    cat("\nGenes of interest:\n")
    print(res_df[res_df$gene %in% genes_of_interest,
                 c("gene", "log2FoldChange", "pvalue", "padj")])
    
    # Save
    csv_file <- paste0(output_prefix, "_", case_label, "_vs_", 
                        control_label, ".csv")
    rds_file <- paste0("results/deseq2/dds_", case_label, "_vs_",
                        control_label, ".rds")
    
    write.csv(res_df, csv_file, row.names = FALSE)
    saveRDS(dds, rds_file)
    
    cat("Saved:", csv_file, "\n")
    cat("Saved:", rds_file, "\n")
    
    return(res_df)
}


# -----------------------------------------------------------------------------
# Identify non-T1 domain mutations (residues 120-216)
# -----------------------------------------------------------------------------
non_t1_samples <- t1_mutations %>%
    filter(residue < 120 | residue > 216)

non_t1_mutant_samples <- non_t1_samples$sample

cat("\nNon-T1 domain mutant samples:\n")
print(non_t1_samples)
cat("Number of non-T1 domain mutant samples:", nrow(non_t1_samples), "\n")

# Quick check - should add up to all mutants
cat("\nSanity check:\n")
cat("T1 mutants:", length(t1_mutant_samples), "\n")
cat("Non-T1 mutants:", length(non_t1_mutant_samples), "\n")
cat("Total mutants:", length(mutant_in_counts), "\n")
cat("T1 + non-T1:", 
    length(t1_mutant_samples) + length(non_t1_mutant_samples), "\n")

# Save
saveRDS(non_t1_mutant_samples, "data/cleaned/non_t1_mutant_samples.rds")

cat("\nNon-T1 domain mutant samples (residues > 216):\n")
print(non_t1_samples)
cat("Number of non-T1 domain mutant samples:", nrow(non_t1_samples), "\n")

# Save
saveRDS(non_t1_mutant_samples, "data/cleaned/non_t1_mutant_samples.rds")


# -----------------------------------------------------------------------------
# Run DESeq2 analyses
# -----------------------------------------------------------------------------

# Analysis 1 — All mutants vs WT
res_all_mutant <- run_deseq2(
    counts_raw     = counts_raw,
    case_samples   = mutant_in_counts,
    control_samples = wt_in_counts,
    case_label     = "allMutant",
    control_label  = "WT",
    output_prefix  = "results/deseq2/DESeq2_KCNA5"
)

# Analysis 2 — T1 mutants vs WT
res_t1_mutant <- run_deseq2(
    counts_raw      = counts_raw,
    case_samples    = t1_mutant_samples,
    control_samples = wt_in_counts,
    case_label      = "T1Mutant",
    control_label   = "WT",
    output_prefix   = "results/deseq2/DESeq2_KCNA5"
)


# -----------------------------------------------------------------------------
# DESeq2 — Non-T1 mutants vs WT
res_non_t1_mutant <- run_deseq2(
    counts_raw      = counts_raw,
    case_samples    = non_t1_mutant_samples,
    control_samples = wt_in_counts,
    case_label      = "nonT1Mutant",
    control_label   = "WT",
    output_prefix   = "results/deseq2/DESeq2_KCNA5"
)
