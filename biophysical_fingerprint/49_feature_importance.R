## ---------------------------------------------------------------------------
## Which of the 24 biophysical dimensions actually drives the cancer-vs-
## background and reference-vs-background separation? A multivariate
## classifier (logistic regression/random forest) would badly overfit at
## n=13 (cancer) or n=35 (reference) against 24 features -- univariate
## standardised mean difference is the safer, more interpretable choice at
## this sample size: for each feature, (group mean - background mean) /
## background sd, ranked by |effect size|. This also lets each finding be
## directly traced to one named biophysical quantity rather than a PCA
## loading (a linear combination across all 24).
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
pc_df <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "full_landscape_pca_umap.csv"))
pc_df$category[pc_df$group == "Observed (COSMIC)"] <- "Observed (COSMIC)"
pc_df$category[pc_df$group == "Background (full saturation)"] <- "Background"
## codon_reachable excluded here -- it's a background-pool membership flag
## (44's dual-null design), not a biophysical feature; leaving it in
## setdiff would silently add a nonsense "feature" row to the effect-size
## table (it slipped through undetected on an earlier run only because
## that run predated this column existing in the saved CSV).
FEATURES <- setdiff(names(pc_df), c("label", "group", "category", "codon_reachable", "PC1", "PC2", "UMAP1", "UMAP2"))

standardised_effect_size <- function(group_vals, bg_vals) {
  (mean(group_vals) - mean(bg_vals)) / sd(bg_vals)
}

compute_effects <- function(group_label, idx_bg) {
  idx_group <- which(pc_df$category == group_label)
  sapply(FEATURES, function(f) standardised_effect_size(pc_df[[f]][idx_group], pc_df[[f]][idx_bg]))
}

idx_bg_full <- which(pc_df$category == "Background")
## The codon-reachable subset -- the statistically meaningful null for
## Observed (COSMIC) specifically (44's dual-null design: cancer mutations
## are real, codon-constrained somatic events, so the biologically relevant
## comparison is against other real mutations that could have occurred, not
## against everything biophysically engineerable). Literature categories
## stay on the full-saturation background, correctly, since they're not
## codon-constrained.
codon_reachable_labels <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "pore_background_substitutions.csv"))$label
idx_bg_codon <- idx_bg_full[pc_df$label[idx_bg_full] %in% codon_reachable_labels]

groups_to_test <- c("C-type_Inactivation", "Activation_Gate", "Drug_Toxin_Binding")
effect_tbl <- do.call(cbind, lapply(groups_to_test, compute_effects, idx_bg = idx_bg_full))
colnames(effect_tbl) <- groups_to_test
cancer_codon <- compute_effects("Observed (COSMIC)", idx_bg_codon)
cancer_full  <- compute_effects("Observed (COSMIC)", idx_bg_full)
effect_df <- as.data.frame(effect_tbl) %>%
  mutate(`Cancer (vs codon-reachable, primary)` = cancer_codon,
         `Cancer (vs full saturation, secondary)` = cancer_full) %>%
  tibble::rownames_to_column("feature") %>%
  select(feature, `Cancer (vs codon-reachable, primary)`, `Cancer (vs full saturation, secondary)`, everything()) %>%
  arrange(desc(abs(`Cancer (vs codon-reachable, primary)`)))

write.csv(effect_df, file.path(RESULTS_DIR, "tables", ANALYSIS, "feature_importance_effect_sizes.csv"), row.names = FALSE)
log_progress("Top 10 features by |standardised effect size| for cancer vs codon-reachable background (primary):")
print(head(effect_df, 10), row.names = FALSE, digits = 2)

## Plot: top 12 features by cancer effect size, grouped bars across the 4 groups
top_features <- head(effect_df$feature, 12)
plot_df <- effect_df %>% filter(feature %in% top_features) %>%
  tidyr::pivot_longer(-feature, names_to = "group", values_to = "effect_size") %>%
  mutate(feature = factor(feature, levels = rev(top_features)))

p <- ggplot(plot_df, aes(x = feature, y = effect_size, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey40") +
  coord_flip() +
  labs(x = NULL, y = "Standardised effect size vs. background\n(group mean - background mean) / background SD",
       title = "Which features distinguish each group from its appropriate background?",
       subtitle = paste0("Top 12 features ranked by |effect size| for cancer vs. codon-reachable background (n=", length(idx_bg_codon), ").\n",
                          "Literature categories and cancer's full-saturation comparison use the full saturation background (n=", length(idx_bg_full), ")."),
       fill = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "feature_importance.png"), p, width = 10, height = 7, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", ANALYSIS, "feature_importance.pdf"), p, width = 10, height = 7)
log_progress("Saved feature_importance_effect_sizes.csv, feature_importance.png/pdf")
