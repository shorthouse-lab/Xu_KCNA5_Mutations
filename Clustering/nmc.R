library(dplyr)
library(tidyr)
library(iPAC)

args <- commandArgs(trailingOnly = TRUE)

gene_name <- "KCNA5"
protein_length <- 613


print(paste0("Running clustering for gene ", gene_name, "..."))
print(paste0("Protein length: ", protein_length, "..."))


# Read file and select for target columns, remove rows with duplicate sample name and mutation_ID
Gene_Data <- read.delim("Clustering/data/cosmic_LI_AC_filtered.tsv") %>%
  filter(MUTATION_DESCRIPTION == "missense_variant" & GENE_SYMBOL == gene_name) %>%
  distinct(SAMPLE_NAME, GENOMIC_MUTATION_ID, .keep_all = TRUE) %>%
  dplyr::select(SAMPLE_NAME, GENOMIC_MUTATION_ID, MUTATION_AA)

# Obtain position of residue by stripping off the wildtype residue and mutant residue from MUTATION_AA
Gene_Data$MUTATION_AA <- substring(Gene_Data$MUTATION_AA, 3)
Gene_Data$MUTATION_POS <- substr(
  Gene_Data$MUTATION_AA,
  2,
  nchar(Gene_Data$MUTATION_AA) - 1
)

samples <- unique(Gene_Data$SAMPLE_NAME)

print("Creating mutation matrix...")

# Create a matrix with numbers of samples as rows and numbers of columns as 
mutation_mat <- matrix(0, nrow = length(samples), ncol = protein_length,
                       dimnames = list(samples, as.character(1:protein_length)))

# Fill the matrix by marking mutation positions with 1
for(i in 1:nrow(Gene_Data)) {
  sample <- Gene_Data$SAMPLE_NAME[i]
  pos <- Gene_Data$MUTATION_POS[i]
  if (as.integer(pos) > protein_length){
    print (c("out of length:", pos))}
  else{
    mutation_mat[sample, pos] <- 1}
}
rownames(mutation_mat) <- NULL

print("Start clustering...")
cluster_result <- nmc(mutation_mat, alpha = 0.05, multtest = "Bonferroni")
write.csv(cluster_result,file=paste0("results/",gene_name,"_clusters.csv"))


