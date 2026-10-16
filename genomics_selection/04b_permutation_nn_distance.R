## ---------------------------------------------------------------------------
## Direct permutation test on the 3D structure: mean nearest-neighbour (NN)
## distance among mutated residues (observed) vs. a null distribution built
## by randomly permuting which residues are "mutated" (same n, same
## background set), matching the original task's own description of
## Analysis 1 ("observed mean nearest-neighbour distance between mutated
## residues versus the permuted null distribution"). Complements
## 04_spacepac_clustering.R's SpaceClust sphere-search approach -- simpler,
## much faster, independent script (doesn't touch the SpacePAC run).
## ---------------------------------------------------------------------------
library(bio3d)
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
NPERM <- 10000
set.seed(42)

cosmic <- readRDS(file.path(PROCESSED_DIR, "structural_cosmic_mutations.rds"))

pdb <- read.pdb("KCNA5.pdb")
ca_sel <- atom.select(pdb, chain = "A", elety = "CA")
ca_atoms <- pdb$atom[ca_sel$atom, ]
coords <- data.frame(resno = ca_atoms$resno, x = ca_atoms$x, y = ca_atoms$y, z = ca_atoms$z)
log_progress("Chain A CA coordinates: ", nrow(coords), " residues (", min(coords$resno), "-", max(coords$resno), ")")

## Mean nearest-neighbour distance among a set of residue numbers
mean_nn_distance <- function(resno_set, coord_df) {
  pts <- coord_df %>% filter(resno %in% resno_set) %>% select(x, y, z) %>% as.matrix()
  if (nrow(pts) < 2) return(NA_real_)
  d <- as.matrix(dist(pts))
  diag(d) <- Inf
  mean(apply(d, 1, min))
}

run_permutation_test <- function(background_resnos, mutated_resnos, coord_df, label) {
  observed <- mean_nn_distance(mutated_resnos, coord_df)
  n_mut <- length(mutated_resnos)
  log_progress(label, ": n_mutated=", n_mut, " of ", length(background_resnos),
               " background residues, observed mean NN distance=", round(observed, 2), " Angstrom")

  null_dist <- replicate(NPERM, {
    perm_set <- sample(background_resnos, n_mut, replace = FALSE)
    mean_nn_distance(perm_set, coord_df)
  })

  ## One-sided: is the observed distance SMALLER (more clustered) than expected by chance?
  p_value <- (sum(null_dist <= observed) + 1) / (NPERM + 1)
  z_score <- (observed - mean(null_dist)) / sd(null_dist)
  log_progress(label, ": null mean=", round(mean(null_dist), 2), " sd=", round(sd(null_dist), 2),
               ", z=", round(z_score, 3), ", p(clustered)=", signif(p_value, 3))

  list(observed = observed, null_dist = null_dist, p_value = p_value, z_score = z_score,
       n_mutated = n_mut, n_background = length(background_resnos), label = label)
}

## --- T1 domain -------------------------------------------------------------
t1_background <- coords$resno[coords$resno >= T1_DOMAIN_START & coords$resno <= T1_DOMAIN_END]
t1_mutated <- unique(cosmic$t1_domain$position)
res_t1 <- run_permutation_test(t1_background, t1_mutated, coords, "T1 domain (128-216)")

## --- Whole modelled protein (120-527) --------------------------------------
whole_background <- coords$resno
whole_mutated <- unique(cosmic$whole_protein$position)
whole_mutated <- whole_mutated[whole_mutated %in% whole_background]  # drop any outside modelled range
res_whole <- run_permutation_test(whole_background, whole_mutated, coords, "Whole modelled protein (120-527)")

## --- Save + plot ------------------------------------------------------------
summary_tab <- bind_rows(lapply(list(res_t1, res_whole), function(r) {
  data.frame(region = r$label, n_mutated = r$n_mutated, n_background = r$n_background,
             observed_mean_nn_dist = r$observed, null_mean = mean(r$null_dist),
             null_sd = sd(r$null_dist), z_score = r$z_score, p_value = r$p_value)
}))
save_table(summary_tab, "permutation_nn_distance_summary", ANALYSIS)
log_progress("Summary:\n", paste(capture.output(print(summary_tab, digits = 3)), collapse = "\n"))

plot_df <- bind_rows(
  data.frame(region = res_t1$label, null_dist = res_t1$null_dist),
  data.frame(region = res_whole$label, null_dist = res_whole$null_dist)
)
obs_df <- data.frame(region = c(res_t1$label, res_whole$label),
                      observed = c(res_t1$observed, res_whole$observed))

p <- ggplot(plot_df, aes(x = null_dist)) +
  geom_histogram(bins = 60, fill = "grey70", colour = NA) +
  geom_vline(data = obs_df, aes(xintercept = observed), colour = "firebrick", linewidth = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  labs(x = "Mean nearest-neighbour distance (Angstrom)", y = "Permutations",
       title = "Observed vs permuted null: mean NN distance among mutated residues",
       subtitle = "Red line = observed value (smaller = more spatially clustered)") +
  theme_kcna5()
save_fig(p, "permutation_nn_distance", ANALYSIS, width = 7, height = 8)

saveRDS(list(res_t1 = res_t1, res_whole = res_whole, coords = coords, summary_tab = summary_tab),
        file.path(PROCESSED_DIR, "structural_permutation_nn_results.rds"))
log_progress("Permutation NN-distance test complete.")
