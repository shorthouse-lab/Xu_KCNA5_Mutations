## ---------------------------------------------------------------------------
## Enumeration for the T1-interface-region FOLDING stability saturation
## screen -- all 19 non-WT substitutions at each position, matching
## 37_pore_full_saturation_enumeration.R's convention for the pore domain.
## This is the missing half of the binding-vs-folding question raised in
## chat: residue 167's contact map (10/11 neighbours same-chain, only 1
## marginal cross-chain contact) reads like an intra-subunit fold-
## stabilising residue rather than an interface one, and its binding ddG is
## indeed near-zero in both orientations -- but we had no FOLDING ddG at all
## for the T1 domain to confirm the other side of that story.
##
## Range extended from the T1_DOMAIN_START/END constants (128-216) down to
## 120-216 -- the interface residue set (interface_residue_positions.csv,
## from the AD/DA binding ddG neighbour searches) spans 120-216, and the
## observed COSMIC mutation R121/121 sits in that 120-127 gap (already
## found and manually patched in chat before this rerun); rather than leave
## the background screen missing exactly the region some observed
## mutations sit in, extending it properly gives a complete, honest
## background for every position that could plausibly matter here.
## ---------------------------------------------------------------------------
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
AA_LIST <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
REGION_START <- 120  # matches interface_residue_positions.csv's own min, not T1_DOMAIN_START (128)
REGION_END <- T1_DOMAIN_END

pdb <- bio3d::read.pdb("data_processed/foldx_mutants/KCNA5_protein_Repair.pdb")
ca_A <- pdb$atom[pdb$atom$chain == "A" & pdb$atom$elety == "CA" &
                    pdb$atom$resno >= REGION_START & pdb$atom$resno <= REGION_END, ]

aa3to1 <- bio3d::aa321(ca_A$resid)
wt_seq <- setNames(aa3to1, ca_A$resno)
log_progress("T1 interface-region WT sequence extracted for ", length(wt_seq), " positions (",
             REGION_START, "-", REGION_END, ")")

enum <- do.call(rbind, lapply(names(wt_seq), function(pos) {
  wt <- wt_seq[[pos]]
  muts <- setdiff(AA_LIST, wt)
  data.frame(position = as.integer(pos), wt_aa = wt, mut_aa = muts,
             label = paste0(wt, pos, muts), stringsAsFactors = FALSE)
}))
enum <- enum[order(enum$position), ]
log_progress("T1 interface-region folding saturation enumeration: ", nrow(enum), " substitutions across ",
             length(wt_seq), " positions")

save_table(enum, "t1_domain_folding_saturation_substitutions", ANALYSIS)
log_progress("Saved t1_domain_folding_saturation_substitutions.csv")
