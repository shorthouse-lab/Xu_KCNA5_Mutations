## ---------------------------------------------------------------------------
## Final consolidated comparison: observed COSMIC cancer mutations (13) vs.
## the full saturation background (1757 -- every codon-reachable-or-not
## substitution possible at these 95 positions, minus the observed/
## reference mutations themselves) vs. the expanded literature reference
## set (35 mutations, 5 mechanistic
## categories: C-type inactivation, selectivity/permeation, activation gate,
## trafficking/biogenesis, drug/toxin binding -- 43_reference_mutation_
## extraction.R), now on 23 dimensions instead of 21: the previous FoldX
## decomposition + HOLE radius + CHAP hydrophobicity + WT structural context,
## plus the two literature-implicated H-bond distances (40_hbond_distances.py):
## D485-W472 (intra-subunit "H-bond timer") and the Y483-T477 inter-subunit
## ring (both instances chain A participates in).
##
## H-bond delta imputation: 75/1805 mutations remove one of these bonds'
## partner atoms entirely (mutating position 472/485/483/477 itself -- e.g.
## D485W, one of our own reference mutations, does exactly this to bond A by
## design). NA in that case doesn't mean "no data", it means "this bond
## cannot exist" -- a more complete disruption than any finite distance
## change, so imputed as 2x the max FINITE |delta| observed for that column
## (a large, unambiguous "broken" sentinel) rather than dropped or set to 0
## (0 would wrongly say "no effect").
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(uwot)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
STRUCT_FEATURES <- c("charge_delta", "volume_delta_A3", "wt_rasa",
                      "wt_flexibility_pctile", "wt_filter_coupling_residual_pctile")
FOLDX_ALL_COLUMNS <- c(
  "ddG_kcal_mol", "backbone_hbond", "sidechain_hbond", "vdw", "electro",
  "solvation_polar", "solvation_hydrophobic", "vdw_clash", "entropy_sidechain",
  "entropy_mainchain", "sloop_entropy", "mloop_entropy", "cis_bond",
  "torsional_clash", "backbone_clash", "helix_dipole", "water_bridge",
  "disulfide", "electrostatic_kon", "partial_covalent", "energy_ionisation",
  "entropy_complex"
)
RADIUS_FEATURE <- "radius_delta_nm"
HYDRO_FEATURE <- "hydro_max_delta"
HBOND_FEATURES <- c("d485_w472_delta", "y483A_t477D_delta", "y483B_t477A_delta")

## --- H-bond distances, all 1805, with NA-broken-bond imputation -----------
hbond_raw <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_hbond_distances.csv"))
hbond_imputed <- hbond_raw
for (col in HBOND_FEATURES) {
  finite_max <- max(abs(hbond_raw[[col]]), na.rm = TRUE)
  n_na <- sum(is.na(hbond_imputed[[col]]))
  hbond_imputed[[col]][is.na(hbond_imputed[[col]])] <- 2 * finite_max * sign(1)  ## always positive-direction "broken" sentinel
  log_progress(col, ": ", n_na, " NA (broken bond) rows imputed as ", round(2 * finite_max, 2))
}
hbond_tbl <- hbond_imputed %>% select(label, all_of(HBOND_FEATURES))

## --- 1. Observed 13 ---------------------------------------------------------
observed <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "mutation_effect_summary.csv")) %>%
  rename(label = mutation) %>%
  select(label, all_of(STRUCT_FEATURES), all_of(RADIUS_FEATURE), all_of(HYDRO_FEATURE)) %>%
  inner_join(read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "observed_foldx_breakdown.csv")), by = "label") %>%
  inner_join(hbond_tbl, by = "label") %>%
  mutate(group = "Observed (COSMIC)")

## --- 2. Background: FULL saturation set (1757), not just codon-reachable --
## Originally this used only the 549 codon-reachable substitutions (the
## "what could a real somatic mutation produce" null) -- correct for
## comparing against the Observed/COSMIC group, since real cancer mutations
## really are codon-constrained. But the Reference (literature) group is
## lab-engineered, NOT codon-constrained (that's exactly why 7 of its 9
## original C-type entries and most of the 35-mutation expansion aren't in
## the codon-reachable set at all -- see chat) -- comparing it against a
## codon-reachable-only null was a mismatch. Switched to the full
## saturation background (every position 433-527 x 19 non-WT AAs, minus the
## observed/reference mutations themselves) so both group comparisons use
## the same, appropriately broad "anything biophysically testable at these
## positions" null. STRUCT_FEATURES (charge/volume/RASA/flexibility/
## filter-coupling) come from pore_full_saturation_substitutions.csv, which
## already covers all 1805 (not just the 1236 that needed fresh FoldX runs).
enum_tbl <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_substitutions.csv")) %>%
  select(label, all_of(STRUCT_FEATURES))
ddg_full <- bind_rows(
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_ddg_progress.csv")) %>% select(label, all_of(FOLDX_ALL_COLUMNS)),
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_ddg_progress.csv")) %>% select(label, all_of(FOLDX_ALL_COLUMNS))
) %>% distinct(label, .keep_all = TRUE)
hole_full <- bind_rows(
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_hole_progress.csv")) %>% select(label, all_of(RADIUS_FEATURE)),
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_hole_progress.csv")) %>% select(label, all_of(RADIUS_FEATURE))
) %>% distinct(label, .keep_all = TRUE)
chap_full <- bind_rows(
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_chap_progress.csv")) %>% select(label, all_of(HYDRO_FEATURE)),
  read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_full_saturation_chap_progress.csv")) %>% select(label, all_of(HYDRO_FEATURE))
) %>% distinct(label, .keep_all = TRUE)

background <- enum_tbl %>%
  inner_join(ddg_full, by = "label") %>%
  inner_join(hole_full, by = "label") %>%
  inner_join(chap_full, by = "label") %>%
  inner_join(hbond_tbl, by = "label") %>%
  filter(!label %in% c(read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "mutation_effect_summary.csv"))$mutation,
                        read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "reference_mutations_extracted.csv"))$label)) %>%
  mutate(group = "Background (full saturation)")
log_progress("Full saturation background: ", nrow(background), " mutations (excludes the 13 observed + 35 reference)")

## --- 3. Reference/literature: 35, 5 mechanistic categories -----------------
reference <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "reference_mutations_extracted.csv")) %>%
  select(label, category, all_of(STRUCT_FEATURES), all_of(FOLDX_ALL_COLUMNS), all_of(RADIUS_FEATURE), all_of(HYDRO_FEATURE)) %>%
  inner_join(hbond_tbl, by = "label") %>%
  mutate(group = "Reference (literature)")

combined <- bind_rows(observed, background, reference)
component_sd <- sapply(combined[setdiff(FOLDX_ALL_COLUMNS, "ddG_kcal_mol")], sd)
FOLDX_COMPONENTS <- names(component_sd[component_sd > 1e-6])
FEATURES <- c(STRUCT_FEATURES, FOLDX_COMPONENTS, RADIUS_FEATURE, HYDRO_FEATURE, HBOND_FEATURES)
stopifnot(!anyNA(combined[FEATURES]))
log_progress("Combined: ", sum(combined$group == "Observed (COSMIC)"), " observed + ",
             sum(combined$group == "Background (full saturation)"), " background + ",
             sum(combined$group == "Reference (literature)"), " reference = ", nrow(combined),
             " total, ", length(FEATURES), " features")

## --- PCA + UMAP -------------------------------------------------------------
X <- scale(as.matrix(combined[FEATURES]))
pca <- prcomp(X, center = FALSE, scale. = FALSE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc_df <- combined %>% mutate(PC1 = pca$x[, 1], PC2 = pca$x[, 2])

set.seed(1)
um <- umap(X, n_neighbors = 15, min_dist = 0.1, metric = "euclidean", n_components = 2)
pc_df <- pc_df %>% mutate(UMAP1 = um[, 1], UMAP2 = um[, 2])

write.csv(pc_df %>% select(label, group, category, PC1, PC2, UMAP1, UMAP2, all_of(FEATURES)),
          file.path(RESULTS_DIR, "tables", ANALYSIS, "full_landscape_pca_umap.csv"), row.names = FALSE)

GROUP_COLOURS <- c("Background (full saturation)" = "grey75",
                    "Observed (COSMIC)" = "firebrick",
                    "Reference (literature)" = "#1b7837")

plot_group <- function(df, xvar, yvar, xlab, ylab, title, subtitle) {
  ggplot(df, aes(.data[[xvar]], .data[[yvar]])) +
    geom_point(data = filter(df, group == "Background (full saturation)"), colour = "grey75", size = 1.4, alpha = 0.6) +
    geom_point(data = filter(df, group != "Background (full saturation)"), aes(colour = group), size = 2.8) +
    ggrepel::geom_text_repel(data = filter(df, group != "Background (full saturation)"),
                              aes(label = label, colour = group), size = 3, fontface = "bold",
                              max.overlaps = 40, show.legend = FALSE) +
    scale_colour_manual(values = GROUP_COLOURS, breaks = c("Observed (COSMIC)", "Reference (literature)")) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle, colour = NULL) +
    theme_bw(base_size = 13) + theme(legend.position = "top")
}

p_pca <- plot_group(pc_df, "PC1", "PC2",
                     paste0("PC1 (", var_explained[1], "%)"), paste0("PC2 (", var_explained[2], "%)"),
                     "Full landscape: COSMIC cancer mutations, literature reference set, and the codon-reachable background",
                     paste0(length(FEATURES), "-D space: FoldX decomposition, HOLE radius, CHAP hydrophobicity, WT context, 2 H-bond distances"))
ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "full_landscape_pca.png"), p_pca, width = 11, height = 8, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "full_landscape_pca.pdf"), p_pca, width = 11, height = 8)

p_umap <- plot_group(pc_df, "UMAP1", "UMAP2", "UMAP1", "UMAP2",
                      "UMAP: full landscape",
                      "Same 23-D space as the PCA above, non-linear embedding")
ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "full_landscape_umap.png"), p_umap, width = 11, height = 8, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "full_landscape_umap.pdf"), p_umap, width = 11, height = 8)

## --- Permutation comparisons: DUAL NULL for the cancer/observed group -----
## Two null models answer genuinely different questions (see chat), not one
## being a "corrected" version of the other:
##   - Codon-reachable background: are cancer mutations atypical relative to
##     other REAL single-nucleotide somatic mutations that could have
##     occurred at these positions? The biologically appropriate null for
##     asking about selection -- matches the trinucleotide-context logic
##     already used for dN/dS and the region-rate tests elsewhere in this
##     pipeline. Only this null applies to Observed (COSMIC), since real
##     cancer mutations really are codon-constrained.
##   - Full saturation background: is this an absolute biophysical outlier,
##     regardless of whether a real mutation could ever produce it? The
##     correct (and only sensible) null for the Reference (literature) group
##     and its categories, since those are lab-engineered substitutions with
##     no codon-reachability constraint at all -- most of the 35 aren't
##     codon-reachable to begin with.
set.seed(1)
NPERM <- 10000
grp <- combined$group
idx_bg_full <- which(grp == "Background (full saturation)")
idx_obs <- which(grp == "Observed (COSMIC)")
idx_ref <- which(grp == "Reference (literature)")

codon_reachable_labels <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_substitutions.csv"))$label
idx_bg_codon <- idx_bg_full[combined$label[idx_bg_full] %in% codon_reachable_labels]
log_progress("Background pools: ", length(idx_bg_full), " full saturation, ", length(idx_bg_codon), " codon-reachable subset")

centroid <- function(idx) colMeans(X[idx, , drop = FALSE])
cdist <- function(a, b) sqrt(sum((a - b)^2))
obs_centroid <- centroid(idx_obs)

perm_vs_background <- function(idx_group, idx_pool) {
  n <- length(idx_group)
  pool_centroid <- centroid(idx_pool)
  observed_dist <- cdist(centroid(idx_group), pool_centroid)
  null <- replicate(NPERM, cdist(colMeans(X[sample(idx_pool, n), , drop = FALSE]), pool_centroid))
  list(dist = observed_dist, null_mean = mean(null), null_sd = sd(null),
       z = (observed_dist - mean(null)) / sd(null), p = mean(null >= observed_dist))
}
perm_vs_target <- function(idx_group, target_centroid, idx_null_pool) {
  n <- length(idx_group)
  observed_dist <- cdist(centroid(idx_group), target_centroid)
  null <- replicate(NPERM, cdist(colMeans(X[sample(idx_null_pool, n), , drop = FALSE]), target_centroid))
  list(dist = observed_dist, null_mean = mean(null), null_sd = sd(null),
       z = (observed_dist - mean(null)) / sd(null), p_closer = mean(null <= observed_dist))
}

res_obs_vs_codon <- perm_vs_background(idx_obs, idx_bg_codon)   ## primary: selection-relevant
res_obs_vs_full  <- perm_vs_background(idx_obs, idx_bg_full)    ## secondary: absolute-outlier check
res_ref_vs_bg    <- perm_vs_background(idx_ref, idx_bg_full)    ## only sensible null for engineered mutations
res_ref_vs_obs   <- perm_vs_target(idx_ref, obs_centroid, idx_bg_full)

log_progress("Cancer vs CODON-REACHABLE background (primary):    z=", round(res_obs_vs_codon$z,2), " p=", round(res_obs_vs_codon$p,4))
log_progress("Cancer vs FULL SATURATION background (secondary):  z=", round(res_obs_vs_full$z,2), " p=", round(res_obs_vs_full$p,4))
log_progress("Reference vs FULL SATURATION background:           z=", round(res_ref_vs_bg$z,2), " p=", round(res_ref_vs_bg$p,4))
log_progress("Reference vs Cancer (closer than chance?):         z=", round(res_ref_vs_obs$z,2), " p_closer=", round(res_ref_vs_obs$p_closer,4))

three_way <- bind_rows(
  data.frame(comparison = "Cancer (observed) vs codon-reachable background (primary, selection-relevant)",
             n_group = length(idx_obs), n_pool = length(idx_bg_codon),
             distance = res_obs_vs_codon$dist, z = res_obs_vs_codon$z, p = res_obs_vs_codon$p),
  data.frame(comparison = "Cancer (observed) vs full saturation background (secondary, absolute-outlier check)",
             n_group = length(idx_obs), n_pool = length(idx_bg_full),
             distance = res_obs_vs_full$dist, z = res_obs_vs_full$z, p = res_obs_vs_full$p),
  data.frame(comparison = "Reference (literature) vs full saturation background",
             n_group = length(idx_ref), n_pool = length(idx_bg_full),
             distance = res_ref_vs_bg$dist, z = res_ref_vs_bg$z, p = res_ref_vs_bg$p),
  data.frame(comparison = "Reference vs Cancer centroid",
             n_group = length(idx_ref), n_pool = length(idx_bg_full),
             distance = res_ref_vs_obs$dist, z = res_ref_vs_obs$z, p = res_ref_vs_obs$p_closer)
)
write.csv(three_way, file.path(RESULTS_DIR, "tables", ANALYSIS, "full_landscape_permutation.csv"), row.names = FALSE)

## Also persist the codon-reachable pool membership as a column on the saved
## PCA/UMAP table (re-saved here with the extra column), so downstream
## scripts (45/46/49) can select the correct null pool per comparison
## without re-deriving it.
codon_flag <- rep(FALSE, nrow(combined))
codon_flag[idx_bg_codon] <- TRUE
pc_df$codon_reachable <- codon_flag
write.csv(pc_df %>% select(label, group, category, codon_reachable, PC1, PC2, UMAP1, UMAP2, all_of(FEATURES)),
          file.path(RESULTS_DIR, "tables", ANALYSIS, "full_landscape_pca_umap.csv"), row.names = FALSE)
log_progress("Re-saved full_landscape_pca_umap.csv with codon_reachable flag")
log_progress("Saved full_landscape_pca_umap.csv, full_landscape_pca.png/pdf, full_landscape_umap.png/pdf, full_landscape_permutation.csv")
