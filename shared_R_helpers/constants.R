## ---------------------------------------------------------------------------
## Shared constants used across all scripts.
## ---------------------------------------------------------------------------

## Cohort = COADREAD (colon + rectal adenocarcinoma), matching the TCGA
## PanCancer Atlas 2018 study on cBioPortal. There is no COAD-only pan-cancer
## atlas study — "coadread_tcga_pan_can_atlas_2018" combines COAD+READ, and
## that combined, curated cohort (clinical incl. MSI, mutation, CNA calls) is
## the source of truth here. GDC is queried with both underlying projects
## (TCGA-COAD, TCGA-READ) so the raw-counts/TPM/methylation pulls match the
## same sample set.
GDC_PROJECTS <- c("TCGA-COAD", "TCGA-READ")
CBIOPORTAL_STUDY <- "coadread_tcga_pan_can_atlas_2018"

KCNA5_SYMBOL   <- "KCNA5"
KCNA5_ENTREZ   <- "3741"
KCNA5_ENSEMBL  <- "ENSG00000130037"

## T1 (S1-S2 linker / tetramerisation) domain boundaries, per task spec
T1_DOMAIN_START <- 128
T1_DOMAIN_END   <- 216

## cBioPortal PanCancer Atlas clinical data carries only continuous MSI
## scores (MSI_SCORE_MANTIS, MSI_SENSOR_SCORE), no pre-binned categorical
## status. MSI-H thresholds below are the published defaults for each method
## (Kautto et al. 2017 for MANTIS; Niu et al. 2014 for MSIsensor). A case is
## called MSI-H if EITHER score crosses its threshold.
MSI_MANTIS_THRESHOLD  <- 0.4
MSI_SENSOR_THRESHOLD  <- 10

## WES capture size assumption for TMB (mutations / Mb). TCGA exome capture
## kits vary (Agilent SureSelect ~38 Mb is the conventional TCGA default);
## documented here as a single overridable constant.
TMB_CAPTURE_SIZE_MB <- 38

RAW_DIR       <- "data_raw"
PROCESSED_DIR <- "data_processed"
RESULTS_DIR   <- "results"
LOG_DIR       <- "logs"

FIG_WIDTH_IN  <- 7
FIG_HEIGHT_IN <- 5
FIG_DPI       <- 300
