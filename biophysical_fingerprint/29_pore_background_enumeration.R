## ---------------------------------------------------------------------------
## Background for the pore/S5/S6 mutation-effect heatmap (24_mutation_effect_
## summary.py): every amino acid change reachable by a SINGLE nucleotide
## substitution from KCNA5's real reference codon, at every residue 433-613
## (the pipeline's existing "S5-P-S6 pore domain" boundary, 12_regional_
## dnds_domains.R), excluding the 13 substitutions actually observed in
## COSMIC. This is the null set the observed mutations get compared against
## -- not "any of the 20 amino acids", which would include changes no point
## mutation could ever produce, but the same codon-reachability constraint
## already used for the region-level dN/dS analysis (11/12_region_dnds*.R).
##
## Only the "free" columns are computed here (no new structure/FoldX/CHAP
## run needed):
##   - charge_delta, volume_delta_A3: sequence only (formulas copied from
##     24_mutation_effect_summary.py for exact consistency with the observed-
##     13 table)
##   - wt_rasa: already computed genome-wide (rasa_all_residues.csv, from
##     Analysis 3/8's DSSP run) -- just needs extracting for the full 433-613
##     range instead of the 13 positions of interest
##   - wt_flexibility_pctile, wt_filter_coupling_residual_pctile: the ANM
##     (bio3d::nma) is already cached (data_processed/enm_modes_cache.rds,
##     built once in 23_enm_scaffold_filter_coupling.R) -- re-running it here
##     for the full range costs nothing beyond loading the cache; extending
##     to more residues doesn't require rebuilding the model
##
## ddG (FoldX BuildModel) is deliberately NOT computed here -- that's the one
## column that needs a new structure per background substitution and is
## handled by the separate orchestration script (30_pore_background_foldx.sh)
## so this enumeration step stays fast and can be re-run freely.
## ---------------------------------------------------------------------------
library(bio3d)
library(dplyr)
source("R/constants.R")
source("R/progress.R")

ANALYSIS <- "structural"
POREDOMAIN_START <- 433
## KCNA5.pdb (the homology model everything downstream is built on) only has
## resolved coordinates for chain A residues 120-527 -- the C-terminal tail
## beyond 527 has no atoms to mutate/score at all, so 613 (the domain's
## sequence-level boundary used in the dN/dS analysis) is not usable here.
## Verified in-session: read.pdb("KCNA5.pdb") chain-A CA range is 120-527
## (408 residues), matching the "408 modelled chain-A residues" background
## already used elsewhere in this pipeline (Analysis 1a/8a permutation
## tests). All 13 observed pore/S5/S6 mutations (max position 525) fall
## inside this range, so nothing observed is lost by the restriction.
POREDOMAIN_END   <- 527

## --- 1. Codon-reachable substitutions per position, 433-613 ----------------
## Uses the whole-gene (T1-split) RefCDS build, gene "KCNA5_nonT1", since the
## domain-specific RefCDS (kcna5_domains_RefCDS.rda) only covers 433-613
## already and would work equally well here, but nonT1 was already verified
## in-session to reproduce KCNA5's real codons (residue numbering
## c(1:127, 217:613), trailing stop codon at position 525 in the concatenated
## sequence dropped) -- reused for exact consistency with that check.
suppressMessages(library(dndscv))
load("data_raw/kcna5_T1_split_RefCDS.rda")
nonT1 <- RefCDS[[which(sapply(RefCDS, function(x) x$gene_name) == "KCNA5_nonT1")]]
cds_seq <- paste(as.character(nonT1$seq_cds), collapse = "")
codons_all <- substring(cds_seq, seq(1, nchar(cds_seq) - 2, 3), seq(3, nchar(cds_seq), 3))
codons_all <- codons_all[1:524]  ## drop trailing stop codon (TGA)
resnum <- c(1:127, 217:613)
stopifnot(length(resnum) == length(codons_all))
names(codons_all) <- resnum

bases <- c("A", "C", "G", "T")
aa1 <- function(cod) toupper(seqinr::translate(strsplit(cod, "")[[1]]))

reachable_substitutions <- function(cod) {
  wt <- aa1(cod)
  out <- character(0)
  for (pos in 1:3) {
    for (b in bases) {
      if (b == substr(cod, pos, pos)) next
      newcod <- cod
      substr(newcod, pos, pos) <- b
      mt <- aa1(newcod)
      if (mt != wt && mt != "*") out <- c(out, mt)
    }
  }
  unique(out)
}

positions <- POREDOMAIN_START:POREDOMAIN_END
rows <- list()
for (r in positions) {
  cod <- codons_all[[as.character(r)]]
  wt_aa <- aa1(cod)
  for (mt_aa in reachable_substitutions(cod)) {
    rows[[length(rows) + 1]] <- data.frame(position = r, wt_aa = wt_aa, mut_aa = mt_aa)
  }
}
background <- bind_rows(rows)
log_progress("Codon-reachable substitutions, ", POREDOMAIN_START, "-", POREDOMAIN_END, ": ",
             nrow(background), " across ", length(positions), " positions")

## --- 2. Exclude the 13 substitutions actually observed in COSMIC -----------
OBSERVED <- data.frame(
  position = c(449, 455, 464, 470, 474, 496, 497, 498, 504, 514, 514, 516, 525),
  mut_aa   = c("P", "T", "C", "T", "E", "M", "D", "L", "R", "I", "L", "I", "W")
)
background <- background %>%
  anti_join(OBSERVED, by = c("position", "mut_aa"))
log_progress("After excluding the 13 observed substitutions: ", nrow(background), " background rows")

## --- 3. Sequence-only columns: charge and side-chain volume change ---------
## Copied verbatim from 24_mutation_effect_summary.py for exact consistency.
FORMAL_CHARGE <- c(D = -1, E = -1, K = 1, R = 1)
VOLUME_A3 <- c(
  A = 88.6, R = 173.4, N = 114.1, D = 111.1, C = 108.5, Q = 143.8,
  E = 138.4, G = 60.1, H = 153.2, I = 166.7, L = 166.7, K = 168.6,
  M = 162.9, F = 189.9, P = 112.7, S = 89.0, T = 116.1, W = 227.8,
  Y = 193.6, V = 140.0
)
charge_of <- function(aa) ifelse(aa %in% names(FORMAL_CHARGE), FORMAL_CHARGE[aa], 0)
background <- background %>%
  mutate(
    charge_delta = unname(charge_of(mut_aa) - charge_of(wt_aa)),
    volume_delta_A3 = unname(VOLUME_A3[mut_aa] - VOLUME_A3[wt_aa]),
    label = paste0(wt_aa, position, mut_aa)
  )

## --- 4. WT RASA, full 433-613 range (already computed genome-wide) --------
rasa_all <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "rasa_all_residues.csv")) %>%
  filter(chain == "A") %>% select(position, wt_rasa = RASA)
background <- background %>% left_join(rasa_all, by = "position")
log_progress("RASA joined; missing: ", sum(is.na(background$wt_rasa)))

## --- 5. WT ANM flexibility + distance-corrected filter-coupling residual,
##        full 433-613 range (cached model, same method as script 23) -------
pdb <- read.pdb("KCNA5.pdb")
ca_pdb <- trim.pdb(pdb, atom.select(pdb, "calpha"))
modes <- readRDS("data_processed/enm_modes_cache.rds")

resno <- ca_pdb$atom$resno
chain <- ca_pdb$atom$chain
flexA <- fluct.nma(modes)[chain == "A"]
resnoA <- resno[chain == "A"]
flex_pctile <- function(r) round(100 * mean(flexA <= flexA[resnoA == r]), 1)

cm <- dccm(modes)
idxA <- which(chain == "A")
filter_idx <- which(chain == "A" & resno %in% 479:484)
corr_to_filter <- sapply(idxA, function(i) mean(cm[i, filter_idx]))
names(corr_to_filter) <- resnoA

xyzA <- as.matrix(ca_pdb$atom[chain == "A", c("x", "y", "z")])
filter_centroid <- colMeans(xyzA[resnoA %in% 479:484, , drop = FALSE])
dist_to_filter <- sqrt(rowSums(sweep(xyzA, 2, filter_centroid)^2))
names(dist_to_filter) <- resnoA

fit <- loess(corr_to_filter ~ dist_to_filter, span = 0.75)
residual <- corr_to_filter - predict(fit, dist_to_filter)
names(residual) <- resnoA
residual_pctile <- function(r) round(100 * mean(residual <= residual[as.character(r)]), 1)

enm_tbl <- data.frame(
  position = resnoA,
  wt_flexibility_pctile = sapply(resnoA, flex_pctile),
  wt_filter_coupling_residual_pctile = sapply(as.character(resnoA), residual_pctile)
) %>% filter(position %in% positions)

background <- background %>% left_join(enm_tbl, by = "position")
log_progress("ENM columns joined; missing flexibility: ", sum(is.na(background$wt_flexibility_pctile)),
             " missing coupling: ", sum(is.na(background$wt_filter_coupling_residual_pctile)))

## --- 6. Save, with placeholder ddG column for the FoldX step to fill ------
background$ddG_kcal_mol <- NA_real_
background$radius_delta_nm <- NA_real_   ## deliberately left NA -- HOLE/CHAP not run for background (scope decision, see chat)
background$hydro_max_delta <- NA_real_   ## same

out_path <- file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_substitutions.csv")
write.csv(background, out_path, row.names = FALSE)
log_progress("Saved ", out_path, " (", nrow(background), " background substitutions, ddG pending)")
