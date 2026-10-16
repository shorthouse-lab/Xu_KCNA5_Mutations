## ---------------------------------------------------------------------------
## Full saturation-mutagenesis enumeration: every one of the 19 non-WT amino
## acids at every one of the 95 structurally-modelled pore/S5/S6 positions
## (433-527, same boundary as 29_pore_background_enumeration.R -- see that
## script's header for why 613 isn't usable, KCNA5.pdb only resolves
## 120-527). 1805 total substitutions (19 x 95), vs. the 562 that are
## reachable by a single real-codon nucleotide change (29's output) -- this
## is the broader "anything engineerable" set requested in chat, not
## constrained to what a real somatic mutation could produce.
##
## Same free (no FoldX/HOLE/CHAP needed) columns as 29: charge_delta,
## volume_delta_A3 (sequence only), wt_rasa/wt_flexibility_pctile/
## wt_filter_coupling_residual_pctile (already computed per-position,
## independent of which mutant AA).
## ---------------------------------------------------------------------------
library(bio3d)
library(dplyr)
source("R/constants.R")
source("R/progress.R")

ANALYSIS <- "structural"
POREDOMAIN_START <- 433
POREDOMAIN_END   <- 527

suppressMessages(library(dndscv))
load("data_raw/kcna5_T1_split_RefCDS.rda")
nonT1 <- RefCDS[[which(sapply(RefCDS, function(x) x$gene_name) == "KCNA5_nonT1")]]
cds_seq <- paste(as.character(nonT1$seq_cds), collapse = "")
codons_all <- substring(cds_seq, seq(1, nchar(cds_seq) - 2, 3), seq(3, nchar(cds_seq), 3))[1:524]
resnum <- c(1:127, 217:613)
names(codons_all) <- resnum
aa1 <- function(cod) toupper(seqinr::translate(strsplit(cod, "")[[1]]))

AAS <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

positions <- POREDOMAIN_START:POREDOMAIN_END
rows <- list()
for (r in positions) {
  wt_aa <- aa1(codons_all[[as.character(r)]])
  for (mt_aa in setdiff(AAS, wt_aa)) {
    rows[[length(rows) + 1]] <- data.frame(position = r, wt_aa = wt_aa, mut_aa = mt_aa)
  }
}
full_sat <- bind_rows(rows) %>%
  mutate(label = paste0(wt_aa, position, mut_aa))
log_progress("Full saturation set: ", nrow(full_sat), " substitutions across ", length(positions), " positions")

## --- Sequence-only columns (same formulas as 29/24) ------------------------
FORMAL_CHARGE <- c(D = -1, E = -1, K = 1, R = 1)
VOLUME_A3 <- c(
  A = 88.6, R = 173.4, N = 114.1, D = 111.1, C = 108.5, Q = 143.8,
  E = 138.4, G = 60.1, H = 153.2, I = 166.7, L = 166.7, K = 168.6,
  M = 162.9, F = 189.9, P = 112.7, S = 89.0, T = 116.1, W = 227.8,
  Y = 193.6, V = 140.0
)
charge_of <- function(aa) ifelse(aa %in% names(FORMAL_CHARGE), FORMAL_CHARGE[aa], 0)
full_sat <- full_sat %>%
  mutate(charge_delta = unname(charge_of(mut_aa) - charge_of(wt_aa)),
         volume_delta_A3 = unname(VOLUME_A3[mut_aa] - VOLUME_A3[wt_aa]))

## --- WT-only position properties (reuse the 29 output directly) -----------
enum_29 <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_substitutions.csv")) %>%
  distinct(position, wt_rasa, wt_flexibility_pctile, wt_filter_coupling_residual_pctile)
full_sat <- full_sat %>% left_join(enum_29, by = "position")
stopifnot(!anyNA(full_sat$wt_rasa))

## --- Flag which rows are already computed (codon-reachable background/
## observed, or the 7 literature ones) so the FoldX runner can skip them ----
already_done <- unique(c(
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_ddg_progress.csv"))$label,
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "mutation_effect_summary.csv"))$mutation,
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "literature_ctype_mutations.csv"))$label
))
full_sat <- full_sat %>% mutate(already_computed = label %in% already_done)
log_progress("Already computed (codon-reachable + observed + literature): ", sum(full_sat$already_computed))
log_progress("Remaining to compute: ", sum(!full_sat$already_computed))

out_path <- file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_substitutions.csv")
write.csv(full_sat, out_path, row.names = FALSE)
log_progress("Saved ", out_path)
