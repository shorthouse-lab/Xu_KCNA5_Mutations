## Analysis 6e (biophysical follow-up): elastic network model (ANM, bio3d::nma)
## on the full WT tetramer, testing two things a static single-mutant
## comparison (6a) can't:
##   (i)  are G497/G504 -- glycine "hinges" per Analysis 8's own framing --
##        predicted to be unusually flexible in the WT network, consistent
##        with a backbone-flexibility role rather than a side-chain one?
##   (ii) are the "stability-only" scaffold residues (S449, F464, A470,
##        A474 -- destabilising per FoldX ddG but geometrically silent per
##        HOLE/CHAP, see Analysis 6b) mechanically *coupled* to the
##        filter's own motion in the WT network, despite not lining it?
##        This is a direct test of Analysis 8's proposed mechanism
##        (scaffold-mediated, indirect effect on filter dynamics/C-type
##        inactivation) that doesn't require moving a single atom via
##        mutation -- it's a property of the WT network topology alone.
##
## Caveat up front: standard ANM/GNM is built purely from the Calpha contact
## topology. A point mutation (side-chain swap) barely moves Calpha
## positions, so comparing WT-vs-mutant normal modes directly would mostly
## just be re-discovering the insensitivity already documented in 6a, for a
## less interesting reason. This analysis instead asks a WT-only question:
## where do these residues sit in the channel's own intrinsic dynamics.

library(bio3d)

RESULTS_DIR <- "results"
ANALYSIS <- "structural"

pdb <- read.pdb("KCNA5.pdb")
ca_sel <- atom.select(pdb, "calpha")
ca_pdb <- trim.pdb(pdb, ca_sel)

cache <- "data_processed/enm_modes_cache.rds"
if (file.exists(cache)) {
  cat("Loading cached ANM modes...\n")
  modes <- readRDS(cache)
} else {
  cat("Building ANM on", nrow(ca_pdb$atom), "CA atoms (4 chains)...\n")
  modes <- nma(ca_pdb, ff = "calpha")
  saveRDS(modes, cache)
}

resno <- ca_pdb$atom$resno
chain <- ca_pdb$atom$chain

## ---------------------------------------------------------------------------
## (i) Per-residue flexibility (chain A), flagging positions of interest
## ---------------------------------------------------------------------------
flex <- fluct.nma(modes)
flexA <- flex[chain == "A"]
resnoA <- resno[chain == "A"]

positions_of_interest <- c(
  S449 = 449, A455 = 455, F464 = 464, A470 = 470, A474 = 474, V496 = 496,
  G497 = 497, S498 = 498, G504 = 504, V514 = 514, V516 = 516, R525 = 525
)
other_glycines <- setdiff(resnoA[which(ca_pdb$atom$resid[chain == "A"] == "GLY")],
                           positions_of_interest)

flex_pctile <- function(r) round(100 * mean(flexA <= flexA[resnoA == r]), 1)

cat("\n--- WT predicted flexibility (ANM), pore-region positions of interest ---\n")
flex_tbl <- data.frame(
  position = names(positions_of_interest),
  resno = positions_of_interest,
  flexibility = sapply(positions_of_interest, function(r) flexA[resnoA == r]),
  flexibility_pctile = sapply(positions_of_interest, flex_pctile)
)
print(flex_tbl, row.names = FALSE)

cat("\n--- for context: all other GLY residues in chain A (n=", length(other_glycines), ") ---\n")
gly_flex <- sapply(other_glycines, function(r) flexA[resnoA == r])
cat("median flexibility percentile of other glycines:",
    round(median(sapply(other_glycines, flex_pctile)), 1), "\n")
cat("G497 percentile:", flex_tbl$flexibility_pctile[flex_tbl$position == "G497"],
    " | G504 percentile:", flex_tbl$flexibility_pctile[flex_tbl$position == "G504"], "\n")

write.csv(flex_tbl, file.path(RESULTS_DIR, "tables", ANALYSIS, "enm_flexibility_pore_positions.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## (ii) Dynamic cross-correlation with the filter (chain A, 479-484)
## ---------------------------------------------------------------------------
cm <- dccm(modes)  # full cross-correlation matrix, all chains x all chains

filter_idx <- which(chain == "A" & resno %in% 479:484)
cat("\nFilter Calpha indices (chain A, 479-484):", length(filter_idx), "atoms\n")

## mean correlation of every chain-A residue with the filter residues
idxA <- which(chain == "A")
corr_to_filter <- sapply(idxA, function(i) mean(cm[i, filter_idx]))
names(corr_to_filter) <- resno[idxA]

corr_pctile <- function(r) round(100 * mean(corr_to_filter <= corr_to_filter[as.character(r)]), 1)

cat("\n--- mean DCCM correlation with filter (chain A 479-484), positions of interest ---\n")
corr_tbl <- data.frame(
  position = names(positions_of_interest),
  resno = positions_of_interest,
  mean_corr_to_filter = sapply(positions_of_interest, function(r) corr_to_filter[as.character(r)]),
  corr_pctile_vs_all_chainA = sapply(positions_of_interest, corr_pctile)
)
print(corr_tbl, row.names = FALSE)

write.csv(corr_tbl, file.path(RESULTS_DIR, "tables", ANALYSIS, "enm_filter_coupling_pore_positions.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## (ii-b) Distance-corrected coupling: raw correlation-to-filter trivially
## rises for any residue simply by being physically closer to the filter in a
## connected elastic network. Fit correlation ~ distance across ALL chain-A
## residues and take residuals for our 11 positions -- a positive residual
## means "more coupled than distance alone predicts", which is the real test
## of Analysis 8's "special scaffold role" claim, not just proximity.
## ---------------------------------------------------------------------------
xyzA <- matrix(ca_pdb$xyz[1, ], ncol = 3, byrow = TRUE)[chain == "A", ]
filter_centroid <- colMeans(xyzA[resnoA %in% 479:484, , drop = FALSE])
dist_to_filter <- sqrt(rowSums(sweep(xyzA, 2, filter_centroid)^2))
names(dist_to_filter) <- resnoA

fit <- loess(corr_to_filter ~ dist_to_filter, span = 0.75)
predicted <- predict(fit, dist_to_filter)
residual <- corr_to_filter - predicted
names(residual) <- resnoA

resid_pctile <- function(r) round(100 * mean(residual <= residual[as.character(r)]), 1)

corr_tbl$distance_to_filter_A <- sapply(positions_of_interest, function(r) dist_to_filter[as.character(r)])
corr_tbl$distance_corrected_residual <- sapply(positions_of_interest, function(r) residual[as.character(r)])
corr_tbl$residual_pctile <- sapply(positions_of_interest, resid_pctile)

cat("\n--- distance-corrected filter coupling (positive = more coupled than distance predicts) ---\n")
print(corr_tbl[, c("position", "resno", "distance_to_filter_A", "mean_corr_to_filter",
                    "distance_corrected_residual", "residual_pctile")], row.names = FALSE)

write.csv(corr_tbl, file.path(RESULTS_DIR, "tables", ANALYSIS, "enm_filter_coupling_pore_positions.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## Figure: flexibility + filter-coupling percentile for each position, split
## by the 6d mechanistic classes (direct/G504R; stability-only; weak-signal)
## ---------------------------------------------------------------------------
library(ggplot2)

class_map <- c(S449 = "stability-only", A455 = "weak signal", F464 = "stability-only",
                A470 = "stability-only", A474 = "stability-only", V496 = "weak signal",
                G497 = "weak signal", S498 = "weak signal", G504 = "direct (G504R)",
                V514 = "weak signal", V516 = "weak signal", R525 = "weak signal")

plot_df <- merge(flex_tbl, corr_tbl, by = c("position", "resno"))
plot_df$class <- class_map[plot_df$position]

## y-axis is the *distance-corrected* residual, not raw correlation -- raw
## correlation with the filter is confounded by simple 3D proximity (every
## position here is already closer to the filter than a random chain-A
## residue), so plotting it directly would visually favour A474/G497 (close,
## high raw correlation) over R525/V516/V514 (far, but genuinely coupled
## *beyond* what their distance predicts) -- the opposite of this section's
## actual finding.
p <- ggplot(plot_df, aes(x = flexibility_pctile, y = residual_pctile,
                          colour = class, label = position)) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey60") +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 4, show.legend = FALSE) +
  labs(x = "WT ANM flexibility percentile (chain A)",
       y = "Distance-corrected filter-coupling residual, percentile (chain A)\n(above 50 = more coupled to filter than distance alone predicts)",
       colour = "Mechanistic class (Analysis 6d)",
       title = "Elastic network model: intrinsic flexibility vs. distance-corrected filter coupling",
       subtitle = "WT tetramer only -- tests whether 'weak signal' residues (6d) are\nmechanically coupled to filter dynamics beyond simple proximity") +
  theme_bw(base_size = 13)

ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "enm_flexibility_filter_coupling.png"),
       p, width = 8, height = 6.5, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "enm_flexibility_filter_coupling.pdf"),
       p, width = 8, height = 6.5)

cat("\nSaved enm_flexibility_pore_positions.csv, enm_filter_coupling_pore_positions.csv,",
    "enm_flexibility_filter_coupling.png/pdf\n")
