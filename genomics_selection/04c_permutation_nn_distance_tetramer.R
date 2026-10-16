## ---------------------------------------------------------------------------
## Tetramer-aware version of 04b: KCNA5 functions as a homotetramer, so a
## somatic mutation alters that residue in ALL FOUR identical subunits at
## once. The single-chain (A-only) NN-distance test in 04b could miss real
## clustering that only appears when a mutated residue's copy on one chain
## sits close to a different mutated residue's copy on an ADJACENT chain,
## at the inter-subunit interface. This redoes the same permutation test
## using all 4 chains' CA coordinates as one combined point cloud.
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
ca_sel <- atom.select(pdb, elety = "CA")
ca_atoms <- pdb$atom[ca_sel$atom, ]
coords <- data.frame(chain = ca_atoms$chain, resno = ca_atoms$resno,
                      x = ca_atoms$x, y = ca_atoms$y, z = ca_atoms$z)
log_progress("All-chain CA coordinates: ", nrow(coords), " atoms across chains ",
             paste(sort(unique(coords$chain)), collapse = ","))

## Mean NN distance for a set of RESIDUE NUMBERS, using ALL chains' copies
## as one combined point cloud (so nearest neighbour can be a different
## chain's copy of the same or a different mutated residue)
mean_nn_distance_tetramer <- function(resno_set, coord_df) {
  pts <- coord_df %>% filter(resno %in% resno_set) %>% select(x, y, z) %>% as.matrix()
  if (nrow(pts) < 2) return(NA_real_)
  d <- as.matrix(dist(pts))
  diag(d) <- Inf
  mean(apply(d, 1, min))
}

run_permutation_test <- function(background_resnos, mutated_resnos, coord_df, label) {
  observed <- mean_nn_distance_tetramer(mutated_resnos, coord_df)
  n_mut <- length(mutated_resnos)
  log_progress(label, ": n_mutated positions=", n_mut, " (x4 chains = ", n_mut * 4, " points), observed mean NN distance=",
               round(observed, 2), " Angstrom")

  null_dist <- replicate(NPERM, {
    perm_set <- sample(background_resnos, n_mut, replace = FALSE)
    mean_nn_distance_tetramer(perm_set, coord_df)
  })

  p_value <- (sum(null_dist <= observed) + 1) / (NPERM + 1)
  z_score <- (observed - mean(null_dist)) / sd(null_dist)
  log_progress(label, ": null mean=", round(mean(null_dist), 2), " sd=", round(sd(null_dist), 2),
               ", z=", round(z_score, 3), ", p(clustered)=", signif(p_value, 3))

  list(observed = observed, null_dist = null_dist, p_value = p_value, z_score = z_score,
       n_mutated = n_mut, label = label)
}

t1_background <- unique(coords$resno[coords$resno >= T1_DOMAIN_START & coords$resno <= T1_DOMAIN_END])
t1_mutated <- unique(cosmic$t1_domain$position)
res_t1 <- run_permutation_test(t1_background, t1_mutated, coords, "T1 domain, all 4 chains (128-216)")

whole_background <- unique(coords$resno)
whole_mutated <- unique(cosmic$whole_protein$position)
whole_mutated <- whole_mutated[whole_mutated %in% whole_background]
res_whole <- run_permutation_test(whole_background, whole_mutated, coords, "Whole protein, all 4 chains (120-527)")

## --- Also report: for each mutated T1 residue, its single closest        ---
## --- INTER-CHAIN approach to any other mutated residue (min across all   ---
## --- chain pairs) -- a more targeted look at interface-spanning contacts ---
interchain_min <- function(mutated_resnos, coord_df) {
  pts <- coord_df %>% filter(resno %in% mutated_resnos)
  d <- as.matrix(dist(pts[, c("x", "y", "z")]))
  diag(d) <- Inf
  same_chain <- outer(pts$chain, pts$chain, "==")
  same_res <- outer(pts$resno, pts$resno, "==")
  d[same_chain] <- Inf  # only keep inter-chain pairs
  data.frame(chain = pts$chain, resno = pts$resno, min_interchain_dist = apply(d, 1, min))
}
interchain_t1 <- interchain_min(t1_mutated, coords) %>% arrange(min_interchain_dist)
log_progress("Closest inter-chain contacts among T1-mutated residues (top 10):\n",
             paste(capture.output(print(head(interchain_t1, 10), row.names = FALSE)), collapse = "\n"))

save_table(bind_rows(
  data.frame(region = res_t1$label, n_mutated = res_t1$n_mutated, observed = res_t1$observed,
             null_mean = mean(res_t1$null_dist), z = res_t1$z_score, p_value = res_t1$p_value),
  data.frame(region = res_whole$label, n_mutated = res_whole$n_mutated, observed = res_whole$observed,
             null_mean = mean(res_whole$null_dist), z = res_whole$z_score, p_value = res_whole$p_value)
), "permutation_nn_distance_tetramer_summary", ANALYSIS)
save_table(interchain_t1, "T1_mutated_closest_interchain_contacts", ANALYSIS)

plot_df <- bind_rows(
  data.frame(region = res_t1$label, null_dist = res_t1$null_dist),
  data.frame(region = res_whole$label, null_dist = res_whole$null_dist)
)
obs_df <- data.frame(region = c(res_t1$label, res_whole$label), observed = c(res_t1$observed, res_whole$observed))
p <- ggplot(plot_df, aes(x = null_dist)) +
  geom_histogram(bins = 60, fill = "grey70", colour = NA) +
  geom_vline(data = obs_df, aes(xintercept = observed), colour = "firebrick", linewidth = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  labs(x = "Mean nearest-neighbour distance (Angstrom, all 4 chains)", y = "Permutations",
       title = "Tetramer-aware: observed vs permuted null NN distance",
       subtitle = "Red line = observed. Nearest neighbour can now be a different chain's copy.") +
  theme_kcna5()
save_fig(p, "permutation_nn_distance_tetramer", ANALYSIS, width = 7, height = 8)

saveRDS(list(res_t1 = res_t1, res_whole = res_whole, interchain_t1 = interchain_t1, coords = coords),
        file.path(PROCESSED_DIR, "structural_permutation_nn_tetramer_results.rds"))
log_progress("Tetramer-aware permutation NN-distance test complete.")
