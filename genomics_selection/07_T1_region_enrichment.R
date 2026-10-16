## ---------------------------------------------------------------------------
## Is the T1 domain, as a whole, enriched for mutations relative to the rest
## of the protein -- at the region level, and at the genuinely 3D/tetramer
## level?
##
## A. Region-level enrichment (Fisher's exact): T1 (128-216) vs rest of the
##    MODELLED protein (120-527, since only the 3D-resolved region is a fair
##    background for a structural analysis). NOTE: because all 4 chains are
##    identical, this ratio is mathematically the same whether computed on
##    one chain or the full tetramer -- the tetramer doesn't change this
##    particular test, and we say so rather than silently re-running a
##    redundant 4x-scaled version.
##
## B. Genuinely tetramer-native 3D test: the four T1 domains fold together
##    to form a compact assembly hub at one end of the channel. Using all 4
##    chains, compute that hub's 3D centroid, then test whether mutated
##    residues (whole-protein) sit closer to it than the background --
##    directly analogous to the existing mutation_simulation.ipynb's
##    selectivity-filter-proximity test, applied to the T1 hub instead.
## ---------------------------------------------------------------------------
library(bio3d)
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
cosmic <- readRDS(file.path(PROCESSED_DIR, "structural_cosmic_mutations.rds"))

pdb <- read.pdb("KCNA5.pdb")
ca_all <- pdb$atom[atom.select(pdb, elety = "CA")$atom, ]
ca_A   <- ca_all %>% filter(chain == "A")

## =============================================================================
## A. Region-level Fisher's exact: T1 vs rest of modelled protein (120-527)
## =============================================================================
log_progress("--- A. T1 domain region-level enrichment (single-chain ratio; ",
             "identical across chains since KCNA5 is a homotetramer) ---")

modelled_positions <- unique(ca_A$resno)
t1_positions <- modelled_positions[modelled_positions >= T1_DOMAIN_START & modelled_positions <= T1_DOMAIN_END]
nonT1_positions <- setdiff(modelled_positions, t1_positions)

whole_mut_positions <- unique(cosmic$whole_protein$position)
whole_mut_positions <- whole_mut_positions[whole_mut_positions %in% modelled_positions]

a <- sum(t1_positions %in% whole_mut_positions)         # T1, mutated
b <- length(t1_positions) - a                             # T1, not mutated
c <- sum(nonT1_positions %in% whole_mut_positions)        # non-T1, mutated
d <- length(nonT1_positions) - c                           # non-T1, not mutated
tab <- matrix(c(a, c, b, d), nrow = 2, dimnames = list(c("T1", "Non-T1"), c("Mutated", "Not mutated")))
log_progress("2x2 table:\n", paste(capture.output(print(tab)), collapse = "\n"))
ft <- fisher.test(tab)
log_progress("T1 vs rest-of-modelled-protein: OR=", round(unname(ft$estimate), 3),
             " 95%CI=[", round(ft$conf.int[1], 3), ", ", round(ft$conf.int[2], 3), "]",
             " p=", signif(ft$p.value, 3),
             " (", a, "/", length(t1_positions), " T1 positions mutated = ",
             round(100 * a / length(t1_positions), 1), "%, vs ", c, "/", length(nonT1_positions),
             " non-T1 = ", round(100 * c / length(nonT1_positions), 1), "%)")

save_table(as.data.frame(tab) %>% tibble::rownames_to_column("region"), "T1_region_enrichment_table", ANALYSIS)
save_table(data.frame(odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
                       p_value = ft$p.value, pct_T1_mutated = 100 * a / length(t1_positions),
                       pct_nonT1_mutated = 100 * c / length(nonT1_positions)),
           "T1_region_enrichment_test", ANALYSIS)

## =============================================================================
## B. Tetramer-native: distance to the T1 assembly's 3D centroid (all 4 chains)
## =============================================================================
log_progress("--- B. Distance to T1 tetramer-assembly centroid (all 4 chains) ---")

ca_t1_allchains <- ca_all %>% filter(resno >= T1_DOMAIN_START, resno <= T1_DOMAIN_END)
t1_centroid <- colMeans(ca_t1_allchains[, c("x", "y", "z")])
log_progress("T1 assembly centroid (all 4 chains, ", nrow(ca_t1_allchains), " atoms): ",
             paste(round(t1_centroid, 1), collapse = ", "))

## Distance from each chain-A residue to the T1 hub centroid
ca_A <- ca_A %>% mutate(
  dist_to_T1_hub = sqrt((x - t1_centroid["x"])^2 + (y - t1_centroid["y"])^2 + (z - t1_centroid["z"])^2)
)
ca_A$is_mutated <- ca_A$resno %in% whole_mut_positions

wt_hub <- wilcox.test(dist_to_T1_hub ~ is_mutated, data = ca_A)
log_progress("Distance-to-T1-hub, mutated vs all modelled residues, Wilcoxon: W=", wt_hub$statistic,
             " p=", signif(wt_hub$p.value, 3),
             " (median mutated=", round(median(ca_A$dist_to_T1_hub[ca_A$is_mutated]), 1),
             "A, median all=", round(median(ca_A$dist_to_T1_hub), 1), "A)")

save_table(ca_A %>% select(resno, resid, dist_to_T1_hub, is_mutated), "T1_hub_distance_all_residues", ANALYSIS)
save_table(data.frame(W = wt_hub$statistic, p_value = wt_hub$p.value), "T1_hub_distance_wilcoxon_test", ANALYSIS)

p_hub <- ggplot(ca_A, aes(x = is_mutated, y = dist_to_T1_hub, fill = is_mutated)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_x_discrete(labels = c("FALSE" = "Not mutated", "TRUE" = "Mutated")) +
  labs(x = NULL, y = "Distance to T1 tetramer-assembly centroid (Angstrom)",
       title = "Are whole-protein mutations closer to the T1 hub than expected?",
       subtitle = paste0("Wilcoxon p=", signif(wt_hub$p.value, 3))) +
  theme_kcna5() + theme(legend.position = "none")
save_fig(p_hub, "T1_hub_distance_by_mutation_status", ANALYSIS)

saveRDS(list(fisher_table = tab, fisher_test = ft, t1_centroid = t1_centroid,
             ca_A = ca_A, wt_hub = wt_hub),
        file.path(PROCESSED_DIR, "structural_T1_region_enrichment_results.rds"))
log_progress("T1 region enrichment analysis complete.")
