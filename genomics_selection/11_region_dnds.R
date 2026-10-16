## ---------------------------------------------------------------------------
## Region-specific dN/dS: is the T1 domain under stronger positive selection
## than the rest of KCNA5, properly correcting for local (trinucleotide-
## context) background mutation rate rather than assuming uniform mutation
## probability across the gene (as the earlier Fisher's-exact region test
## did)?
##
## Built via dndscv::buildref(), splitting KCNA5 into two synthetic "genes"
## (KCNA5_T1, KCNA5_nonT1) sharing the real CDS but partitioned by residue
## range -- an explicitly supported use of buildref() per its own vignette
## ("you could annotate different protein domains under a different gene
## name to perform inferences of selection at a domain level").
##
## GENOME BUILD: our COSMIC data is GRCh37/hg19 (verified empirically --
## KCNA5 mutation genomic coordinates match Ensembl's GRCh37 gene coordinates
## exactly, not GRCh38, despite the Xu repo's cosmic_dnds_calculation.R
## specifying refdb="hg38"). Using hg19 throughout, which is also dndscv's
## own default/bundled build.
##
## KCNA5 CDS (verified via Ensembl + cross-checked against 3 real COSMIC
## mutation coordinates, all landing in the correct codon):
##   chr12:5,153,314-5,155,155, GRCh37, + strand, phase 0, single exon,
##   1842nt = 613 aa + stop. Residue N starts at 5153314 + (N-1)*3.
## ---------------------------------------------------------------------------
library(dndscv)
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
## Real genomic coordinates (GRCh37), used only to map mutation positions
CDS_START_GENOMIC <- 5153314
CDS_END_GENOMIC   <- 5155155
## /tmp/kcna5_grch37_cds.fa only contains the small excerpt chr12:5153313-5155156
## (not the whole chromosome), so buildref() must be given coordinates LOCAL to
## that fasta (1-based from its first base), not real chromosome coordinates --
## buildref validates cds coordinates against the fasta's sequence length, and
## real genomic coordinates (~5.15 million) exceed that length, which is
## exactly the "coordinates fall outside chromosome length" error this fixes.
FASTA_GENOMIC_START <- 5153313
to_local <- function(genomic_pos) genomic_pos - FASTA_GENOMIC_START + 1
CDS_START <- to_local(CDS_START_GENOMIC)
CDS_END   <- to_local(CDS_END_GENOMIC)

## =============================================================================
## A. Build the custom cdsfile (T1 vs non-T1 pseudo-genes)
## =============================================================================
res_to_gpos <- function(res) CDS_START + (res - 1) * 3  # local pos of codon start

t1_g_start <- res_to_gpos(T1_DOMAIN_START)
t1_g_end   <- res_to_gpos(T1_DOMAIN_END) + 2
nonT1_seg1_end   <- t1_g_start - 1                 # residues 1..(T1_START-1)
nonT1_seg2_start <- t1_g_end + 1                   # residues (T1_END+1)..613+stop, to CDS_END
nonT1_len1 <- nonT1_seg1_end - CDS_START + 1
nonT1_len2 <- CDS_END - nonT1_seg2_start + 1
nonT1_total_len <- nonT1_len1 + nonT1_len2

log_progress("T1 CDS segment: ", t1_g_start, "-", t1_g_end, " (", t1_g_end - t1_g_start + 1, " nt)")
log_progress("non-T1 segment 1 (residues 1-", T1_DOMAIN_START - 1, "): ", CDS_START, "-", nonT1_seg1_end,
             " (", nonT1_len1, " nt)")
log_progress("non-T1 segment 2 (residues ", T1_DOMAIN_END + 1, "-end): ", nonT1_seg2_start, "-", CDS_END,
             " (", nonT1_len2, " nt)")
stopifnot((t1_g_end - t1_g_start + 1) == (T1_DOMAIN_END - T1_DOMAIN_START + 1) * 3)
stopifnot((nonT1_total_len + (t1_g_end - t1_g_start + 1)) == (CDS_END - CDS_START + 1))

## The two KCNA5_nonT1 rows are the two exon-like segments of ONE CDS (the
## gene split around the T1 domain), so they must share a single cds.id --
## dndscv's buildref groups rows by cds.id to reassemble one coding
## sequence (as multi-exon genes like PIK3CA do in dndscv's own bundled
## example data). Giving them different cds.id values instead makes
## buildref treat them as two competing alternative transcripts and keep
## only one (verified empirically: doing so silently dropped every
## non-T1 mutation from the analysis, in the earlier failed run).
cdsfile_df <- data.frame(
  gene.id = c("KCNA5_T1", "KCNA5_nonT1", "KCNA5_nonT1"),
  gene.name = c("KCNA5_T1", "KCNA5_nonT1", "KCNA5_nonT1"),
  cds.id = c("KCNA5_T1_cds1", "KCNA5_nonT1_cds1", "KCNA5_nonT1_cds1"),
  chr = c("12", "12", "12"),
  chr.coding.start = c(t1_g_start, CDS_START, nonT1_seg2_start),
  chr.coding.end   = c(t1_g_end, nonT1_seg1_end, CDS_END),
  cds.start = c(1, 1, nonT1_len1 + 1),
  cds.end   = c(t1_g_end - t1_g_start + 1, nonT1_len1, nonT1_total_len),
  length = c(t1_g_end - t1_g_start + 1, nonT1_total_len, nonT1_total_len),
  strand = c(1, 1, 1)
)
log_progress("Custom cdsfile:\n", paste(capture.output(print(cdsfile_df)), collapse = "\n"))

cdsfile_path <- file.path(RAW_DIR, "kcna5_T1_split_cdsfile.txt")
write.table(cdsfile_df, cdsfile_path, sep = "\t", row.names = FALSE, quote = FALSE)

## =============================================================================
## B. Build the custom RefCDS via buildref()
## =============================================================================
genome_fasta <- "/tmp/kcna5_grch37_cds.fa"  # chr12:5153313-5155156, GRCh37 (1bp margin each side for trinucleotide context)
refcds_out <- file.path(RAW_DIR, "kcna5_T1_split_RefCDS.rda")

log_progress("Running buildref()...")
buildref(cdsfile = cdsfile_path, genomefile = genome_fasta, outfile = refcds_out, useids = FALSE)
log_progress("Custom RefCDS built: ", refcds_out)

## =============================================================================
## C. Prepare mutation data: ALL KCNA5 mutation types (not just missense --
##    dndscv needs synonymous mutations to calibrate the local neutral rate)
## =============================================================================
log_progress("--- Preparing full KCNA5 mutation set (all types, GRCh37) ---")
raw <- read.csv("/tmp/kcna5_cosmic_raw.csv", stringsAsFactors = FALSE)
log_progress("Raw KCNA5 COSMIC records: ", nrow(raw))
log_progress("Mutation types:\n", paste(capture.output(print(table(raw$MUTATION_DESCRIPTION))), collapse = "\n"))

## dndscv wants sampleID, chr, pos, ref, mut -- one row per (sample, variant)
## pos must be converted to the same LOCAL coordinate system as the custom
## RefCDS (built against the small fasta excerpt, not real chromosome
## coordinates) -- otherwise every mutation position would fall far outside
## the gene models and dndscv would find nothing coding.
mut_data <- raw %>%
  filter(MUTATION_DESCRIPTION %in% c("missense_variant", "synonymous_variant", "stop_gained", "stop_lost")) %>%
  distinct(SAMPLE_NAME, CHROMOSOME, GENOME_START, GENOMIC_MUT_ALLELE, .keep_all = TRUE) %>%
  transmute(sampleID = SAMPLE_NAME, chr = as.character(CHROMOSOME), pos = to_local(as.integer(GENOME_START)),
            ref = GENOMIC_WT_ALLELE, mut = GENOMIC_MUT_ALLELE)
log_progress("Mutations for dndscv (missense+synonymous+nonsense, deduped): ", nrow(mut_data))
stopifnot(all(mut_data$pos >= CDS_START & mut_data$pos <= CDS_END))
save_table(mut_data, "kcna5_dnds_input_mutations", ANALYSIS)

## =============================================================================
## D. Run dndscv with the custom T1/non-T1 reference
## =============================================================================
## outp=3 (default) would additionally fit dNdScv's negative-binomial
## regression across genes to model rate overdispersion using covariates --
## a step designed for exome-wide analyses with thousands of genes, and
## which fails outright (theta.ml, "missing value where TRUE/FALSE needed")
## with only 2 genes to regress over. outp=2 stops at dNdSloc: a per-gene
## Poisson MLE of dN/dS with its own likelihood-based CI, requiring no
## cross-gene regression -- exactly the level appropriate for comparing
## two hand-defined pseudo-genes rather than genome-wide gene discovery.
log_progress("--- Running dndscv() with custom RefCDS (outp=2: dNdSloc, no cross-gene NB regression) ---")
dndsout <- dndscv(mut_data, refdb = refcds_out, cv = NULL, outp = 2,
                   max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

sel_cv <- dndsout$sel_loc
log_progress("Region-level dN/dS results (dNdSloc, per-gene MLE):\n", paste(capture.output(print(sel_cv)), collapse = "\n"))

## Approximate 95% CI for wmis_loc via the delta method on the log rate
## ratio. dndscv computes wmis_loc = (n_mis/exp_mis) / (n_syn/exp_syn)
## (its "mrfold" term is exactly n_syn/exp_syn) -- exp_mis/exp_syn are
## fixed, trinucleotide-model-derived quantities, so essentially all the
## sampling variance is in the two Poisson counts n_mis and n_syn, giving
## Var(log wmis_loc) ~= 1/n_mis + 1/n_syn (standard log-rate-ratio CI).
sel_cv <- sel_cv %>%
  mutate(
    se_log_wmis = sqrt(1 / n_mis + 1 / n_syn),
    wmis_loc_lo = wmis_loc * exp(-1.96 * se_log_wmis),
    wmis_loc_hi = wmis_loc * exp( 1.96 * se_log_wmis)
  )
log_progress("With approximate 95% CIs on wmis_loc:\n",
             paste(capture.output(print(sel_cv[, c("gene_name","n_syn","n_mis","wmis_loc","wmis_loc_lo","wmis_loc_hi","pmis_loc","qmis_loc")], digits = 3)), collapse = "\n"))
save_table(sel_cv, "kcna5_T1_vs_nonT1_dnds", ANALYSIS)

## =============================================================================
## E. Plot: missense dN/dS with CI, T1 vs non-T1
## =============================================================================
library(ggplot2)
sel_cv$gene_name <- factor(sel_cv$gene_name, levels = c("KCNA5_T1", "KCNA5_nonT1"))
p <- ggplot(sel_cv, aes(x = gene_name, y = wmis_loc)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(ymin = wmis_loc_lo, ymax = wmis_loc_hi), size = 0.8, colour = "firebrick") +
  geom_text(aes(label = paste0("n_mis=", n_mis, ", n_syn=", n_syn, "\np=", signif(pmis_loc, 2))),
            vjust = -1.2, size = 3.2) +
  labs(x = NULL, y = "Missense dN/dS (dNdSloc MLE, 95% CI)",
       title = "Region-specific dN/dS: T1 domain vs rest of KCNA5",
       subtitle = "dndscv::buildref() custom two-gene split, trinucleotide-context corrected background") +
  theme_kcna5()
save_fig(p, "region_dnds_T1_vs_nonT1", ANALYSIS, width = 6, height = 5.5)

saveRDS(list(dndsout = dndsout, cdsfile_df = cdsfile_df, mut_data = mut_data, sel_cv = sel_cv),
        file.path(PROCESSED_DIR, "structural_region_dnds_results.rds"))
log_progress("Region-specific dN/dS analysis complete.")
