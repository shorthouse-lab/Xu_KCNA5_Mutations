## ---------------------------------------------------------------------------
## Pore-radius profiling (HOLE-style, via CHAP -- the Channel Annotation
## Package, Klesse et al. 2019) to look at where mutations sit relative to
## the actual conduction pathway, particularly the pore-domain (S5-P-S6)
## mutations flagged by the region-level dN/dS analysis (12).
##
## CHAP isn't packaged for this machine's toolchain (old, 2020, linux-64-only
## bioconda build) so it's run inside a small Docker image built for this
## analysis (kcna5-chap:0.9.1, from mambaorg/micromamba + `conda install
## bioconda::chap`, run under linux/amd64 emulation since the host is
## arm64) -- see scratchpad chap_docker/Dockerfile. Rebuild with:
##   docker build --platform linux/amd64 -t kcna5-chap:0.9.1 <dockerfile_dir>
##
## KEY METHODOLOGICAL NOTES (verified empirically before trusting any of
## this, since none of it can be taken for granted with a homology model
## that was never built with any particular orientation in mind):
##
## 1. Pore axis: KCNA5.pdb is a raw homology model with no guarantee its
##    coordinate frame aligns to any membrane normal. Confirmed by fitting
##    a plane through the 4 chains' CA centroids (a homotetramer's 4
##    subunits are arranged with their centroids in a plane perpendicular
##    to the pore/symmetry axis) that the true axis runs along global
##    +y (0.0134, 0.9999, 0.0088), NOT CHAP's default (0,0,1) -- residual
##    of the plane fit is ~0.12 Angstrom, confirming a very clean 4-fold
##    symmetric arrangement and a reliable axis estimate. Passed explicitly
##    via -pf-chan-dir-vec; CHAP's default z-axis assumption would have
##    (and initially did) send the probe sideways through the protein core.
##
## 2. Units: CHAP/GROMACS work internally in nm; the PDB (and this script's
##    other structural scripts) are in Angstrom. All coordinates passed to
##    -pf-init-probe-pos must be divided by 10 -- an easy silent error
##    (gives "pore radius at initial probe position is infinite", since the
##    probe ends up placed ~10x too far from the protein to find any atoms).
##
## 3. Validation: the profile's true minimum-radius location was checked
##    two independent ways before trusting it: (a) it lands at residue
##    ~484, immediately downstream of the sequence-identified canonical Kv
##    selectivity-filter motif T477-T478-V479-G480-Y481-G482(-D483), a
##    hallmark K+-channel signature found by direct regex search of chain
##    A's sequence, entirely independent of the structural calculation;
##    (b) the model's own pre-placed HETATM K ions (4 of them, evenly
##    spaced ~2.5-3.5 Angstrom apart along the same axis) sit in the same
##    general region. The *magnitude* of the minimum radius found here
##    (~0.19-0.26nm) is noticeably wider than typical experimental HOLE
##    profiles of K+ channel selectivity filters (~0.10-0.15nm) -- almost
##    certainly because homology modelling (side-chain packing from a
##    template) does not reliably reproduce the precise backbone geometry
##    that makes real selectivity filters so tight. Treat radius values
##    near the filter as directionally informative, not literally precise.
##
## 4. Residue <-> pore-position mapping: CHAP's actual path-finding method
##    (inplane_optim, the HOLE algorithm) optimises the probe's position
##    *within* each perpendicular plane to maximise local radius, so the
##    true centreline can bend off the straight init-pos/direction-vector
##    line assumed here. Rather than try to recover that bent centreline,
##    every residue's own CA atom is projected onto the straight pore axis
##    to get its axial coordinate s_residue, and the profile's radius(s) is
##    read off at that same s. This asks a well-posed question ("how tight
##    is the channel at the height where this residue sits") without
##    requiring the exact 3D path shape, and is the standard way pore-radius
##    profiles get overlaid with sequence/domain annotation in the
##    literature.
## ---------------------------------------------------------------------------
library(bio3d)
library(dplyr)
library(ggplot2)
library(jsonlite)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
CHAP_IMAGE <- "kcna5-chap:0.9.1"
INIT_POS_NM <- c(13.2852, -13.28204, 13.31009)   # tetramer CA centroid, nm
DIR_VEC     <- c(0.01344207, 0.9998711, 0.008781797)  # pore axis (plane-normal of 4 chain centroids)

chap_dir <- file.path(PROCESSED_DIR, "chap")
dir.create(chap_dir, showWarnings = FALSE)

## =============================================================================
## A. Run CHAP (skip if already run)
## =============================================================================
chap_json <- file.path(chap_dir, "kcna5_pore.json")
if (!file.exists(chap_json)) {
  log_progress("Running CHAP pore-profile analysis via Docker...")
  wd <- normalizePath(chap_dir)
  cmd <- sprintf(
    'docker run --rm --platform linux/amd64 -v %s:/data -v %s:/out -w /out %s -f /data/KCNA5.pdb -s /data/KCNA5.pdb -sel-pathway Protein -pf-chan-dir-vec %s -pf-init-probe-pos %s -pf-max-free-dist 5 -out-filename kcna5_pore',
    shQuote(normalizePath(".")), shQuote(wd), CHAP_IMAGE,
    paste(DIR_VEC, collapse = " "), paste(INIT_POS_NM, collapse = " ")
  )
  system(cmd)
  stopifnot(file.exists(chap_json))
}
log_progress("CHAP output: ", chap_json)

## =============================================================================
## B. Load profile, build residue -> axial-position -> local-radius mapping
##
## SUPERSEDED, 2026-07: this straight-axis-projection lookup into CHAP's own
## radiusMin(s) profile was found to disagree sharply -- and specifically --
## with CHAP's own OBJ mesh (both its geometry and its colour coding) and
## with a from-scratch nearest-atom check, at several axial heights (found
## via R525, s~1.385nm: this method said 0.076nm, everything else said
## ~0.85-0.9nm; also affected residue 234 and the whole pore/selectivity
## filter domain's aggregate ranking). Structure mismatch was directly ruled
## out (identical atom count, identical R525 coordinates, all CHAP outputs
## from one docker run within an 8-second window). Best explanation: CHAP's
## in-plane radius optimiser (a stochastic simulated-annealing search per
## its own documented parameters) failing to converge at specific heights.
## chap_residue_pore_radius_map.csv and chap_mutated_residues_pore_radius.csv
## are now authoritatively produced by 19_direct_pore_radius.py instead,
## which recomputes local radius directly from atomic coordinates (widest
## inscribed sphere via lateral grid search, same axis/probe parameters
## validated below, no dependency on CHAP's own radius search). Re-running
## this section will silently regress the two tables back to the flawed
## values -- run 19_direct_pore_radius.py afterward if this script is rerun.
## =============================================================================
chap_out <- fromJSON(chap_json)
prof <- chap_out$pathwayProfile
s_grid <- prof$s
r_grid <- prof$radiusMin

pdb <- read.pdb("KCNA5.pdb")
ca <- pdb$atom[atom.select(pdb, elety = "CA")$atom, ]
ca$s_axial <- as.numeric((as.matrix(ca[, c("x","y","z")]) / 10 - matrix(INIT_POS_NM, nrow(ca), 3, byrow = TRUE)) %*% DIR_VEC)
ca$pore_radius_nm <- approx(s_grid, r_grid, xout = ca$s_axial, rule = 1)$y

## Sanity check against the sequence-validated selectivity filter
sf_check <- ca %>% filter(chain == "A", resno %in% 477:484) %>% select(resno, resid, s_axial, pore_radius_nm)
log_progress("Selectivity filter motif (477-484) sanity check:\n", paste(capture.output(print(sf_check)), collapse = "\n"))

## =============================================================================
## C. Merge with domain boundaries and COSMIC mutation data
## =============================================================================
domains <- tribble(
  ~domain, ~start, ~end,
  "Disordered (N-term)", 1, 120,
  "Linker", 121, 127,
  "T1", 128, 216,
  "Linker", 217, 230,
  "S0", 231, 240,
  "Linker", 241, 246,
  "S1", 247, 270,
  "Disordered", 271, 321,
  "S2", 322, 344,
  "S3", 345, 376,
  "Linker", 377, 391,
  "S4", 392, 416,
  "Linker", 417, 432,
  "S5", 433, 456,
  "Linker", 457, 467,
  "Pore/selectivity filter", 468, 510,
  "S6/C-term", 511, 613
)
classify_domain <- function(pos) {
  hit <- domains[pos >= domains$start & pos <= domains$end, ]
  if (nrow(hit) == 0) return(NA_character_)
  hit$domain[1]
}

## use chain A as the representative subunit (all 4 are ~identical by symmetry)
res_map <- ca %>% filter(chain == "A") %>%
  transmute(resno, resid, s_axial, pore_radius_nm, domain = sapply(resno, classify_domain))
save_table(res_map, "chap_residue_pore_radius_map", ANALYSIS)

cosmic <- readRDS(file.path(PROCESSED_DIR, "structural_cosmic_mutations.rds"))
mut_positions <- unique(cosmic$whole_protein$position)
res_map$is_mutated <- res_map$resno %in% mut_positions

mutated_map <- res_map %>% filter(is_mutated) %>% arrange(pore_radius_nm)
log_progress("Mutated residues ranked by local pore radius (narrowest first):\n",
             paste(capture.output(print(mutated_map)), collapse = "\n"))
save_table(mutated_map, "chap_mutated_residues_pore_radius", ANALYSIS)

## Narrowness percentile: what fraction of ALL modelled residues are narrower
## than this one (so "how close to a constriction" is comparable across
## residues without depending on the model's absolute radius miscalibration)
res_map <- res_map %>% mutate(narrowness_pctile = 100 * rank(pore_radius_nm) / sum(!is.na(pore_radius_nm)))
mutated_map <- res_map %>% filter(is_mutated) %>% arrange(narrowness_pctile)

## =============================================================================
## D. Plot: pore radius profile with domain bands + mutations overlaid
## =============================================================================
domain_bands <- domains %>% filter(!domain %in% c("Linker")) %>%
  rowwise() %>%
  mutate(s_start = min(res_map$s_axial[res_map$resno >= start & res_map$resno <= end], Inf, na.rm = TRUE),
         s_end   = max(res_map$s_axial[res_map$resno >= start & res_map$resno <= end], -Inf, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.finite(s_start), is.finite(s_end))

profile_df <- data.frame(s = s_grid, radius = r_grid) %>% filter(radius < 3)  # drop huge open-vestibule tails

p <- ggplot() +
  geom_rect(data = domain_bands, aes(xmin = s_start, xmax = s_end, ymin = -0.05, ymax = -0.01, fill = domain)) +
  geom_line(data = profile_df, aes(x = s, y = radius), colour = "grey30") +
  geom_hline(yintercept = 0.15, linetype = "dotted", colour = "firebrick") +
  geom_point(data = res_map %>% filter(is_mutated), aes(x = s_axial, y = pore_radius_nm), colour = "firebrick", size = 2) +
  annotate("text", x = -Inf, y = 0.16, label = "~K+ ionic radius (0.15nm)", hjust = -0.05, size = 3, colour = "firebrick") +
  coord_cartesian(xlim = range(res_map$s_axial, na.rm = TRUE) + c(-0.3, 0.3)) +
  labs(x = "Axial position along pore (nm)", y = "Pore radius (nm)", fill = "Domain",
       title = "KCNA5 pore radius profile (CHAP) with mutations overlaid",
       subtitle = "Red points = mutated residues, plotted at their own axial position and CHAP's local pore radius there") +
  theme_kcna5() + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
save_fig(p, "chap_pore_radius_profile", ANALYSIS, width = 11, height = 6)

## =============================================================================
## E. Focused plot: mutated residues, ranked by narrowness percentile
## =============================================================================
mutated_map_plot <- mutated_map %>% mutate(resid_label = paste0(resid, resno))
p2 <- ggplot(mutated_map_plot, aes(x = reorder(resid_label, narrowness_pctile), y = narrowness_pctile, fill = domain)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Narrowness percentile (0 = at the tightest constriction, 100 = widest)",
       title = "Mutated KCNA5 residues ranked by proximity to the pore's tightest point",
       fill = "Domain") +
  theme_kcna5()
save_fig(p2, "chap_mutated_residues_narrowness_rank", ANALYSIS, width = 8, height = 10)

saveRDS(list(chap_out = chap_out, res_map = res_map, mutated_map = mutated_map),
        file.path(PROCESSED_DIR, "structural_chap_pore_results.rds"))
log_progress("CHAP pore-radius analysis complete.")
