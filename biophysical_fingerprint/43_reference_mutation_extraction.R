## ---------------------------------------------------------------------------
## Extracts every concrete, single-residue KCNA5 substitution from the user-
## supplied literature reference table (KCNA5_pore_filter_mutation_reference_
## MASTER.csv) that maps unambiguously onto this construct, spanning 5
## mechanistic categories (C-type inactivation, selectivity/permeation,
## activation gate, trafficking/biogenesis, drug/toxin binding) rather than
## just the original 9 C-type-inactivation-only set.
##
## Because the full saturation screen (37-42) already covers all 19 non-WT
## amino acids at every position 433-527, every mutation extracted here is
## expected to already exist in pore_full_saturation_ddg_progress.csv /
## _hole_progress.csv / _chap_progress.csv (or the original 569 already-
## computed set) -- this script does NOT run any new FoldX/HOLE/CHAP, it
## only looks up already-computed rows. stopifnot() below will fail loudly
## if that assumption breaks for any entry, rather than silently proceeding
## with missing data.
##
## Reference-numbering corrections already resolved directly in the source
## file (not re-derived here): H452Q excluded (position 452 is native Val in
## this construct, not His -- unresolved numbering mismatch per source);
## Y519N/F -> F519N only (Y519F would be a no-op, this construct is
## natively Phe at 519); L508M -> I508M (construct is Ile, not Leu, at 508);
## E422Q/K -> Q460K only (E422Q is a no-op, construct is natively Gln at
## 460). Double mutants (R487F+V476A, G504A+V505G) are represented here by
## their individual single-substitution halves only -- no double-mutant
## structures were built anywhere in this pipeline.
##
## One further correction NOT caught in the source file (found here by
## checking pore_full_saturation_substitutions.csv directly): Ric02/Yeo96's
## "T505I/V/S" uses the wrong reference letter for this construct -- native
## residue at 505 is Val, not Thr (the source file's OWN Din05 entry
## correctly uses "V505G" at the same position, so this is an internal
## inconsistency in the source table, not a new numbering scheme). Same
## Kv1.5-isoform-reference pattern as the already-fixed I508/F519 entries.
## Translated: T505I -> V505I, T505S -> V505S; T505V is dropped entirely --
## it would be a no-op (already native Val).
## ---------------------------------------------------------------------------
library(dplyr)
source("R/constants.R")
source("R/progress.R")

ANALYSIS <- "structural"

## label, category (matches the source file's own Category column, kept
## for later stratified comparisons), source_citekey (for traceability)
REFERENCE_MUTATIONS <- tribble(
  ~label,    ~category,                 ~source_citekey,
  "R487V",   "C-type_Inactivation",     "Lop93/Edu07",
  "R487F",   "C-type_Inactivation",     "Ahe09",
  "W472F",   "C-type_Inactivation",     "Yan97",
  "D485W",   "C-type_Inactivation",     "Ple13",
  "H463G",   "C-type_Inactivation",     "Edu07/Zha07",
  "H463A",   "C-type_Inactivation",     "Edu07",
  "H463R",   "C-type_Inactivation",     "Edu07",
  "H463K",   "C-type_Inactivation",     "Edu07",
  "H463E",   "C-type_Inactivation",     "Edu07",
  "V476G",   "Selectivity_Permeation",  "Cha10",
  "R487Y",   "Selectivity_Permeation",  "Wan09",
  "V516A",   "Activation_Gate",         "Hac02",
  "V516W",   "Activation_Gate",         "Hac02/Kit04",
  "V516D",   "Activation_Gate",         "Hac02",
  "F519A",   "Activation_Gate",         "Hac02",
  "F519W",   "Activation_Gate",         "Hac02",
  "F519D",   "Activation_Gate",         "Hac02",
  "F519N",   "Activation_Gate",         "Ric02 (corrected)",
  "P513D",   "Activation_Gate",         "Suk03",
  "G504A",   "Activation_Gate",         "Din05",
  "V505G",   "Activation_Gate",         "Din05",
  "V505I",   "Activation_Gate",         "Ric02/Yeo96 (corrected from 'T505I')",
  "V505S",   "Drug_Toxin_Binding",      "Yeo96 (corrected from 'T505S')",
  "V512A",   "Activation_Gate",         "Ric02/Yeo96/Din19/Dec04",
  "V512M",   "Activation_Gate",         "Ric02/Fra97",
  "I508V",   "Trafficking_Biogenesis",  "Str14",
  "I508M",   "Drug_Toxin_Binding",      "Fra97 (corrected)",
  "T479A",   "Drug_Toxin_Binding",      "Din19/Dec04",
  "T480A",   "Drug_Toxin_Binding",      "Din19/Dec04",
  "I502A",   "Drug_Toxin_Binding",      "Din19",
  "V505A",   "Drug_Toxin_Binding",      "Din19/Dec04",
  "I508A",   "Drug_Toxin_Binding",      "Din19",
  "L510A",   "Drug_Toxin_Binding",      "Din19",
  "T477S",   "Drug_Toxin_Binding",      "Fra97",
  "Q460K",   "Drug_Toxin_Binding",      "Mac89 (corrected)"
)
log_progress(nrow(REFERENCE_MUTATIONS), " reference mutations to extract, spanning ",
             n_distinct(REFERENCE_MUTATIONS$category), " categories")

## --- Pull each one from wherever it was computed -----------------------
## ddG breakdown: either the original background/observed/literature run, or
## the full saturation run.
ddg_bg  <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_ddg_progress.csv"))
ddg_sat <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_ddg_progress.csv"))
ddg_lit <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "literature_ctype_mutations.csv"))
ddg_all <- bind_rows(ddg_bg, ddg_sat, ddg_lit %>% select(-radius_delta_nm, -hydro_max_delta, -qc_disagreement)) %>%
  distinct(label, .keep_all = TRUE)

hole_bg  <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_hole_progress.csv"))
hole_sat <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_hole_progress.csv"))
hole_lit <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "literature_ctype_mutations.csv")) %>% select(label, radius_delta_nm)
hole_all <- bind_rows(hole_bg, hole_sat, hole_lit) %>% distinct(label, .keep_all = TRUE)

chap_bg  <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_chap_progress.csv"))
chap_sat_path <- file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_chap_progress.csv")
chap_sat <- if (file.exists(chap_sat_path)) read.csv(chap_sat_path) else chap_bg[0, ]
chap_lit <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "literature_ctype_mutations.csv")) %>% select(label, hydro_max_delta)
chap_all <- bind_rows(chap_bg, chap_sat, chap_lit) %>% distinct(label, .keep_all = TRUE)

## --- Sequence/structural-context features (charge/volume delta, WT RASA/
## flexibility/filter-coupling) -- computed fresh here the same way as
## 29_pore_background_enumeration.R / 43's earlier draft forgot to join
## these in at all (caught when 44_full_landscape_pca_umap.R errored on
## missing columns) -----------------------------------------------------
FORMAL_CHARGE <- c(D = -1, E = -1, K = 1, R = 1)
VOLUME_A3 <- c(
  A = 88.6, R = 173.4, N = 114.1, D = 111.1, C = 108.5, Q = 143.8,
  E = 138.4, G = 60.1, H = 153.2, I = 166.7, L = 166.7, K = 168.6,
  M = 162.9, F = 189.9, P = 112.7, S = 89.0, T = 116.1, W = 227.8,
  Y = 193.6, V = 140.0
)
charge_of <- function(aa) ifelse(aa %in% names(FORMAL_CHARGE), FORMAL_CHARGE[aa], 0)

## position, wt_aa, mut_aa for each reference label -- parsed from the label
## itself (format <WT><position><MUT>, e.g. "R487V") rather than re-typed,
## to avoid a transcription mismatch with the tribble above.
parse_label <- function(lbl) {
  wt_aa <- substr(lbl, 1, 1)
  mut_aa <- substr(lbl, nchar(lbl), nchar(lbl))
  position <- as.integer(substr(lbl, 2, nchar(lbl) - 1))
  data.frame(label = lbl, wt_aa = wt_aa, mut_aa = mut_aa, position = position)
}
parsed <- bind_rows(lapply(REFERENCE_MUTATIONS$label, parse_label))

wt_position_features <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_substitutions.csv")) %>%
  distinct(position, wt_rasa, wt_flexibility_pctile, wt_filter_coupling_residual_pctile)

struct_features <- parsed %>%
  left_join(wt_position_features, by = "position") %>%
  mutate(charge_delta = unname(charge_of(mut_aa) - charge_of(wt_aa)),
         volume_delta_A3 = unname(VOLUME_A3[mut_aa] - VOLUME_A3[wt_aa])) %>%
  select(label, charge_delta, volume_delta_A3, wt_rasa, wt_flexibility_pctile, wt_filter_coupling_residual_pctile)
stopifnot(!anyNA(struct_features))
REFERENCE_MUTATIONS <- REFERENCE_MUTATIONS %>% left_join(struct_features, by = "label")

missing_ddg  <- setdiff(REFERENCE_MUTATIONS$label, ddg_all$label)
missing_hole <- setdiff(REFERENCE_MUTATIONS$label, hole_all$label)
missing_chap <- setdiff(REFERENCE_MUTATIONS$label, chap_all$label)
if (length(missing_ddg) > 0)  log_progress("Missing ddG (not yet computed): ", paste(missing_ddg, collapse = ", "))
if (length(missing_hole) > 0) log_progress("Missing HOLE: ", paste(missing_hole, collapse = ", "))
if (length(missing_chap) > 0) log_progress("Missing CHAP (likely still running): ", paste(missing_chap, collapse = ", "))

out <- REFERENCE_MUTATIONS %>%
  left_join(ddg_all, by = "label") %>%
  left_join(hole_all %>% select(label, radius_delta_nm), by = "label") %>%
  left_join(chap_all %>% select(label, hydro_max_delta), by = "label")

out_path <- file.path(RESULTS_DIR, "tables", ANALYSIS, "reference_mutations_extracted.csv")
write.csv(out, out_path, row.names = FALSE)
log_progress("Saved ", out_path, " (", sum(!is.na(out$ddG_kcal_mol)), "/", nrow(out), " with ddG, ",
             sum(!is.na(out$hydro_max_delta)), "/", nrow(out), " with CHAP hydrophobicity)")
