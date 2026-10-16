## ---------------------------------------------------------------------------
## Follow-up to 13_chap_pore_analysis.R: the pore-domain (S5-P-S6) mutations
## flagged as under region-level positive selection (12_regional_dnds_domains.R)
## turned out NOT to sit at the pore's actual narrowest constriction (13) --
## they cluster around the 75-80th narrowness percentile, i.e. moderately
## open. If they're not doing something sterically obvious, a natural next
## question is whether they perturb the local pore-lining/pore-facing
## HYDROPHOBICITY environment instead -- CHAP computes this by default
## alongside the radius profile (Klesse et al. 2019's headline extension over
## classic HOLE is exactly this: combining geometry with a physicochemical
## profile to flag putative hydrophobic gates).
##
## Two distinct things are computed here, kept clearly separate:
##  (A) The MODEL's own local pore-facing/lining hydrophobicity at each
##      mutated position (CHAP's pfHydrophobicityMean/plHydrophobicityMean,
##      an average over whichever residues from all 4 chains CHAP judged to
##      be pathway-facing/lining at that axial height, using the *current*
##      WT model) -- describes what kind of local chemical environment a
##      given axial position already has.
##  (B) The DIRECTION AND MAGNITUDE of hydrophobicity change caused by each
##      *specific observed substitution*, computed directly from the exact
##      per-residue Wimley-White 1996 scale values CHAP itself uses
##      internally (pulled straight from the installed
##      chap/share/data/hydrophobicity/wimley_white_1996.json rather than
##      retyped from memory, to guarantee it matches (A) exactly: +1 =
##      most hydrophobic, -1 = most hydrophilic).
##
## Caveat carried over from 13: there is no explicit solvent in this static
## homology model, so CHAP's water-density/dewetting fields are empty here
## -- this analysis is about pore-lining side-chain chemistry, not a direct
## test of hydrophobic gating (which strictly requires MD with explicit
## water to observe wet/dry bistability).
## ---------------------------------------------------------------------------
library(bio3d)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(jsonlite)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"

## =============================================================================
## A. Load CHAP profile + the exact WT/mutant hydrophobicity scale it used
## =============================================================================
chap_out <- fromJSON(file.path(PROCESSED_DIR, "chap", "kcna5_pore.json"))
prof <- chap_out$pathwayProfile
s_grid <- prof$s

ww <- fromJSON("data_raw/chap_scales/wimley_white_1996.json")$hydrophobicity
ww_scale <- setNames(ww$hydrophobicity, ww$resname)
aa3to1 <- c(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G",
            HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S",
            THR="T", TRP="W", TYR="Y", VAL="V")
aa1to3 <- setNames(names(aa3to1), aa3to1)

## =============================================================================
## B. Residue -> axial position -> local model hydrophobicity (reuse the
##    same axis/units established and validated in 13_chap_pore_analysis.R)
## =============================================================================
INIT_POS_NM <- c(13.2852, -13.28204, 13.31009)
DIR_VEC     <- c(0.01344207, 0.9998711, 0.008781797)

pdb <- read.pdb("KCNA5.pdb")
ca <- pdb$atom[atom.select(pdb, elety = "CA")$atom, ]
ca$s_axial <- as.numeric((as.matrix(ca[, c("x","y","z")]) / 10 - matrix(INIT_POS_NM, nrow(ca), 3, byrow = TRUE)) %*% DIR_VEC)
ca$pore_radius_nm <- approx(s_grid, prof$radiusMin, xout = ca$s_axial, rule = 1)$y
ca$pf_hydrophobicity <- approx(s_grid, prof$pfHydrophobicityMean, xout = ca$s_axial, rule = 1)$y
ca$pl_hydrophobicity <- approx(s_grid, prof$plHydrophobicityMean, xout = ca$s_axial, rule = 1)$y

res_map <- ca %>% filter(chain == "A") %>%
  transmute(resno, resid, s_axial, pore_radius_nm, pf_hydrophobicity, pl_hydrophobicity)

## =============================================================================
## C. Pore-domain mutations flagged as under selection (12): actual WT->MUT
##    substitutions, with Delta-hydrophobicity from the exact CHAP scale
## =============================================================================
raw <- read.csv("/tmp/kcna5_cosmic_raw.csv", stringsAsFactors = FALSE)
missense <- raw %>% filter(MUTATION_DESCRIPTION == "missense_variant")
missense$pos <- as.integer(gsub("[^0-9]", "", missense$MUTATION_AA))
missense$wt_aa1  <- substr(missense$MUTATION_AA, 3, 3)
missense$mut_aa1 <- substr(missense$MUTATION_AA, nchar(missense$MUTATION_AA), nchar(missense$MUTATION_AA))

POREDOM_POSITIONS <- c(470, 474, 496, 497, 498, 504)   # from region-level dN/dS follow-up (script 12)

pore_muts <- missense %>%
  filter(pos %in% POREDOM_POSITIONS) %>%
  distinct(SAMPLE_NAME, pos, wt_aa1, mut_aa1, MUTATION_AA) %>%
  mutate(
    wt_hydrophob  = ww_scale[aa1to3[wt_aa1]],
    mut_hydrophob = ww_scale[aa1to3[mut_aa1]],
    delta_hydrophob = mut_hydrophob - wt_hydrophob
  ) %>%
  left_join(res_map, by = c("pos" = "resno")) %>%
  arrange(pos)

log_progress("Pore-domain mutations under selection: WT/mutant hydrophobicity and local pore context:\n",
             paste(capture.output(print(pore_muts %>% select(pos, MUTATION_AA, wt_hydrophob, mut_hydrophob,
                                                               delta_hydrophob, pore_radius_nm, pf_hydrophobicity))), collapse = "\n"))
save_table(pore_muts, "chap_pore_domain_hydrophobicity_changes", ANALYSIS)

## =============================================================================
## D. Plot: Delta-hydrophobicity of each substitution vs local model context
## =============================================================================
pore_muts_u <- pore_muts %>% distinct(pos, MUTATION_AA, .keep_all = TRUE) %>%
  mutate(label = paste0(MUTATION_AA, " (", resid, pos, ")"))

p1 <- ggplot(pore_muts_u, aes(x = reorder(label, pos), y = delta_hydrophob, fill = delta_hydrophob > 0)) +
  geom_col() +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                     labels = c("TRUE" = "more hydrophobic", "FALSE" = "more hydrophilic"), name = NULL) +
  labs(x = NULL, y = "Delta hydrophobicity\n(mutant - WT, Wimley-White 1996, CHAP scale [-1,1])",
       title = "Hydrophobicity change at pore-domain mutations under selection",
       subtitle = "Positive = substitution makes the position more hydrophobic") +
  theme_kcna5() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p1, "chap_pore_domain_delta_hydrophobicity", ANALYSIS, width = 9, height = 5.5)

## Context plot: where do these positions sit on the whole-channel
## hydrophobicity profile (not just their own Delta)?
profile_df <- data.frame(s = s_grid, pf_hydrophob = prof$pfHydrophobicityMean, radius = prof$radiusMin) %>%
  filter(radius < 3)
p2 <- ggplot() +
  geom_line(data = profile_df, aes(x = s, y = pf_hydrophob), colour = "grey30") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(data = pore_muts_u, aes(x = s_axial, y = pf_hydrophobicity), colour = "firebrick", size = 3) +
  geom_text_repel(data = pore_muts_u, aes(x = s_axial, y = pf_hydrophobicity, label = paste0(resid, pos)),
                   size = 3.2, min.segment.length = 0) +
  coord_cartesian(xlim = range(res_map$s_axial[res_map$resno %in% 430:613], na.rm = TRUE)) +
  labs(x = "Axial position along pore (nm)", y = "Pore-facing hydrophobicity (CHAP, [-1,1])",
       title = "Pore-facing hydrophobicity profile around the pore domain",
       subtitle = "Red points = axial position of the 6 pore-domain mutations under selection") +
  theme_kcna5()
save_fig(p2, "chap_pore_domain_hydrophobicity_context", ANALYSIS, width = 9, height = 5.5)

saveRDS(list(pore_muts = pore_muts, res_map = res_map, ww_scale = ww_scale),
        file.path(PROCESSED_DIR, "structural_chap_hydrophobicity_results.rds"))
log_progress("Pore-domain hydrophobicity analysis complete.")
