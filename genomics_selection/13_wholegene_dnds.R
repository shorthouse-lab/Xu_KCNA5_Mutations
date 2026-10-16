## ---------------------------------------------------------------------------
## Whole-gene KCNA5 dN/dS (no domain split) -- companion to 11_region_dnds.R
## and 12_regional_dnds_domains.R, which only ever computed SPLIT dN/dS.
## Same buildref()/dndscv() methodology, GRCh37 coordinates, COSMIC input as
## those scripts, just a single gene model instead of 2 or 4 pseudo-genes.
##
## The original /tmp/kcna5_grch37_cds.fa (1bp-margin CDS excerpt) no longer
## exists on disk (ephemeral /tmp file from an earlier session) -- rebuilt
## here from the CACHED, already-verified RefCDS sequence fragments in
## data_raw/kcna5_T1_split_RefCDS.rda (KCNA5_T1 + KCNA5_nonT1's seq_cds,
## reassembled in genomic order, plus the 1bp up/downstream margin bases
## recovered from nonT1's seq_cds1up/seq_cds1down shifted-sequence fields)
## rather than re-deriving the sequence from any external source. Verified
## before use: reassembled CDS is exactly 1842nt, starts ATG, ends TGA.
## ---------------------------------------------------------------------------
library(dndscv)
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
FASTA_GENOMIC_START <- 5153313
to_local <- function(genomic_pos) genomic_pos - FASTA_GENOMIC_START + 1
CDS_START <- to_local(5153314)  # = 2
CDS_END   <- to_local(5155155)  # = 1843

## =============================================================================
## A. Single-gene cdsfile + buildref()
## =============================================================================
cdsfile_df <- data.frame(
  gene.id = "KCNA5", gene.name = "KCNA5", cds.id = "KCNA5_cds1",
  chr = "12", chr.coding.start = CDS_START, chr.coding.end = CDS_END,
  cds.start = 1, cds.end = CDS_END - CDS_START + 1, length = CDS_END - CDS_START + 1,
  strand = 1
)
cdsfile_path <- file.path(RAW_DIR, "kcna5_wholegene_cdsfile.txt")
write.table(cdsfile_df, cdsfile_path, sep = "\t", row.names = FALSE, quote = FALSE)

genome_fasta <- "/tmp/kcna5_wholegene_cds.fa"
refcds_out <- file.path(RAW_DIR, "kcna5_wholegene_RefCDS.rda")
log_progress("Running buildref() for whole-gene KCNA5...")
buildref(cdsfile = cdsfile_path, genomefile = genome_fasta, outfile = refcds_out, useids = FALSE)

## =============================================================================
## B. Mutation data -- reuse the exact saved input from 11_region_dnds.R
## (same filter: missense/synonymous/stop_gained/stop_lost, deduped, already
## in the correct local-fasta coordinate system)
## =============================================================================
mut_data <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "kcna5_dnds_input_mutations.csv"))
log_progress("Mutations for whole-gene dndscv: ", nrow(mut_data))

## =============================================================================
## C. Run dndscv (dNdSloc, outp=2 -- see 11_region_dnds.R header for why)
## =============================================================================
dndsout <- dndscv(mut_data, refdb = refcds_out, cv = NULL, outp = 2,
                   max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
sel_cv <- dndsout$sel_loc
log_progress("Whole-gene dN/dS result:\n", paste(capture.output(print(sel_cv)), collapse = "\n"))

sel_cv <- sel_cv %>%
  mutate(
    se_log_wmis = sqrt(1 / n_mis + 1 / n_syn),
    wmis_loc_lo = wmis_loc * exp(-1.96 * se_log_wmis),
    wmis_loc_hi = wmis_loc * exp( 1.96 * se_log_wmis)
  )
log_progress("With approximate 95% CI:\n",
             paste(capture.output(print(sel_cv[, c("gene_name","n_syn","n_mis","wmis_loc","wmis_loc_lo","wmis_loc_hi","pmis_loc","qmis_loc")], digits = 3)), collapse = "\n"))

frac <- max(0, (sel_cv$wmis_loc - 1) / sel_cv$wmis_loc)
log_progress("Estimated excess (functionally relevant) fraction of missense mutations: ",
             round(100 * frac, 1), "% (", round(sel_cv$n_mis * frac, 1), " of ", sel_cv$n_mis, ")")

save_table(sel_cv, "kcna5_wholegene_dnds", ANALYSIS)
log_progress("Whole-gene dN/dS analysis complete.")
