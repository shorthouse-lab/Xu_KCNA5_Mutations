## ---------------------------------------------------------------------------
## Region-specific dN/dS across KCNA5's actual functional domains, not just
## a T1-vs-"rest of gene" binary split (11_region_dnds.R). A T1-vs-everything
## comparison buries whatever is happening in the voltage-sensing (S1-S4) and
## pore (S5-P-S6) modules inside one pooled bucket -- the standard way this
## kind of analysis is run is domain-by-domain: T1 (tetramerization), S0,
## the voltage-sensor domain (S1-S4), and the pore domain (S5-P-S6).
##
## Domain boundaries match results/tables/structural/domain_level_mutation_rate.csv
## (built in 08_summary_plots.R), collapsed into the 4 standard functional
## modules rather than each individual TM helix, since several helices have
## too few mutations individually to fit dNdSloc meaningfully:
##   T1                 128-216  (tetramerization domain)
##   S0                  231-240  (S0 helix)
##   VSD (S1-S4)         247-416  (voltage-sensing domain, incl. S1-S2/S2-S3/S3-S4 loops)
##   Pore domain (S5-P-S6) 433-613  (incl. S5-P loop, selectivity filter, S6/C-term)
## Disordered N-terminus and inter-domain linkers (1-127, 217-230, 241-246,
## 417-432, 457-467) are deliberately excluded from all 4 gene models --
## mutations landing there simply aren't assigned to any region (this is
## fine: buildref doesn't require full CDS coverage, and dNdSloc calibrates
## each region's neutral rate from its own local synonymous mutations, not
## a genome/gene-wide pool).
##
## Same GRCh37 CDS, coordinate system, and buildref()/dndscv() methodology
## as 11_region_dnds.R -- see that script for the full genome-build
## verification and buildref mechanics (local-fasta coordinates, shared
## cds.id for multi-segment genes, trailing-newline fix, etc).
## ---------------------------------------------------------------------------
library(dndscv)
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
CDS_START_GENOMIC <- 5153314
CDS_END_GENOMIC   <- 5155155
FASTA_GENOMIC_START <- 5153313
to_local <- function(genomic_pos) genomic_pos - FASTA_GENOMIC_START + 1
CDS_START <- to_local(CDS_START_GENOMIC)
res_to_gpos <- function(res) CDS_START + (res - 1) * 3

## =============================================================================
## A. Define the 4 functional-domain pseudo-genes
## =============================================================================
region_defs <- tribble(
  ~gene.name,             ~res_start, ~res_end,
  "KCNA5_T1",             128,        216,
  "KCNA5_S0",             231,        240,
  "KCNA5_VSD_S1S4",       247,        416,
  "KCNA5_Pore_S5PS6",     433,        613
)

cdsfile_df <- region_defs %>%
  mutate(
    gene.id = gene.name,
    cds.id = paste0(gene.name, "_cds1"),
    chr = "12",
    chr.coding.start = res_to_gpos(res_start),
    chr.coding.end = res_to_gpos(res_end) + 2,
    length = chr.coding.end - chr.coding.start + 1,
    cds.start = 1,
    cds.end = length,
    strand = 1
  ) %>%
  select(gene.id, gene.name, cds.id, chr, chr.coding.start, chr.coding.end, cds.start, cds.end, length, strand)
log_progress("Functional-domain cdsfile:\n", paste(capture.output(print(cdsfile_df)), collapse = "\n"))
stopifnot(all(cdsfile_df$length %% 3 == 0))

cdsfile_path <- file.path(RAW_DIR, "kcna5_domains_cdsfile.txt")
write.table(cdsfile_df, cdsfile_path, sep = "\t", row.names = FALSE, quote = FALSE)

## =============================================================================
## B. Build the custom RefCDS
## =============================================================================
genome_fasta <- "/tmp/kcna5_grch37_cds.fa"
refcds_out <- file.path(RAW_DIR, "kcna5_domains_RefCDS.rda")

log_progress("Running buildref()...")
buildref(cdsfile = cdsfile_path, genomefile = genome_fasta, outfile = refcds_out, useids = FALSE)
log_progress("Custom domain-level RefCDS built: ", refcds_out)

## =============================================================================
## C. Mutation data (identical extraction to 11_region_dnds.R: all coding
##    variant types, deduped, positions converted to the same local coord
##    system as the fasta used to build both RefCDS objects)
## =============================================================================
raw <- read.csv("/tmp/kcna5_cosmic_raw.csv", stringsAsFactors = FALSE)
mut_data <- raw %>%
  filter(MUTATION_DESCRIPTION %in% c("missense_variant", "synonymous_variant", "stop_gained", "stop_lost")) %>%
  distinct(SAMPLE_NAME, CHROMOSOME, GENOME_START, GENOMIC_MUT_ALLELE, .keep_all = TRUE) %>%
  transmute(sampleID = SAMPLE_NAME, chr = as.character(CHROMOSOME), pos = to_local(as.integer(GENOME_START)),
            ref = GENOMIC_WT_ALLELE, mut = GENOMIC_MUT_ALLELE)
log_progress("Mutations for dndscv: ", nrow(mut_data))

## =============================================================================
## D. Run dndscv (dNdSloc, outp=2 -- see 11_region_dnds.R for why outp=3's
##    cross-gene NB regression isn't usable with only a handful of genes)
## =============================================================================
log_progress("--- Running dndscv() across 4 functional domains ---")
dndsout <- dndscv(mut_data, refdb = refcds_out, cv = NULL, outp = 2,
                   max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

sel_loc <- dndsout$sel_loc %>%
  mutate(
    se_log_wmis = sqrt(1 / n_mis + 1 / n_syn),
    wmis_loc_lo = wmis_loc * exp(-1.96 * se_log_wmis),
    wmis_loc_hi = wmis_loc * exp( 1.96 * se_log_wmis)
  ) %>%
  arrange(match(gene_name, c("KCNA5_T1", "KCNA5_S0", "KCNA5_VSD_S1S4", "KCNA5_Pore_S5PS6")))
log_progress("Domain-level dN/dS:\n",
             paste(capture.output(print(sel_loc[, c("gene_name","n_syn","n_mis","wmis_loc","wmis_loc_lo","wmis_loc_hi","pmis_loc","qmis_loc")], digits = 3)), collapse = "\n"))
save_table(sel_loc, "kcna5_domain_dnds", ANALYSIS)

## =============================================================================
## E. Plot
## =============================================================================
label_map <- c(KCNA5_T1 = "T1\n(128-216)", KCNA5_S0 = "S0\n(231-240)",
                KCNA5_VSD_S1S4 = "S1-S4\n(247-416)", KCNA5_Pore_S5PS6 = "S5-P-S6\n(433-613)")
sel_loc$region_label <- factor(label_map[sel_loc$gene_name], levels = label_map)

p <- ggplot(sel_loc, aes(x = region_label, y = wmis_loc)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(ymin = wmis_loc_lo, ymax = wmis_loc_hi), size = 0.8, colour = "firebrick") +
  geom_text(aes(label = paste0("n_mis=", n_mis, ", n_syn=", n_syn, "\np=", signif(pmis_loc, 2))),
            vjust = -1.1, size = 3) +
  labs(x = NULL, y = "Missense dN/dS (dNdSloc MLE, 95% CI)",
       title = "Region-specific dN/dS across KCNA5 functional domains",
       subtitle = "dndscv::buildref() 4-way domain split, trinucleotide-context corrected background") +
  theme_kcna5()
save_fig(p, "region_dnds_by_domain", ANALYSIS, width = 8, height = 5.5)

saveRDS(list(dndsout = dndsout, cdsfile_df = cdsfile_df, mut_data = mut_data, sel_loc = sel_loc),
        file.path(PROCESSED_DIR, "structural_domain_dnds_results.rds"))
log_progress("Domain-level dN/dS analysis complete.")
