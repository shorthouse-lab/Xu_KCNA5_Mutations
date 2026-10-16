## ---------------------------------------------------------------------------
## xCell cross-check of run_immune_infiltration_any_alteration.R's quanTIseq
## result (immune_infiltration_allmutant_vs_wt_quantiseq.csv, AllMutant vs
## WT, MSI+TMB-adjusted linear model -- null throughout, best p_adj=0.80 for
## NK cell). xCell uses a different approach (ssGSEA-based enrichment
## scores against a much larger, non-compositional signature set, 39 cell
## types vs quanTIseq's ~10 fractions that sum to 1) -- running it
## independently asks whether that null result is a genuine absence of
## immune-compositional difference, or specific to quanTIseq's constrained-
## regression assumptions.
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis4_network_immune"
deconv <- readRDS(file.path(PROCESSED_DIR, "analysis4_immunedeconv_results.rds"))
groups <- readRDS(file.path(PROCESSED_DIR, "patient_groups.rds"))
clin <- readRDS(file.path(RAW_DIR, "clinical_cbioportal.rds"))

xc <- deconv$xcell
xc_long <- xc %>%
  ## xCell includes 3 composite/summary scores alongside the 36 cell-type
  ## scores ("immune score", "stroma score", "microenvironment score") --
  ## these are not cell fractions, drop them here so the per-cell-type model
  ## below isn't testing a score built from a mix of other rows in the same
  ## table.
  filter(!cell_type %in% c("immune score", "stroma score", "microenvironment score")) %>%
  tidyr::pivot_longer(-cell_type, names_to = "sample", values_to = "score") %>%
  mutate(patient_id = gsub("\\.", "-", sample)) %>%
  select(-sample)

msi_lookup <- clin %>% transmute(patient_id = patientId, msi_score = as.numeric(MSI_SENSOR_SCORE))

model_df <- xc_long %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = score) %>%
  left_join(groups %>% select(patient_id, group, TMB), by = "patient_id") %>%
  left_join(msi_lookup, by = "patient_id") %>%
  filter(group %in% c("T1_mutant", "Non_T1_mutant", "WT"), !is.na(msi_score), !is.na(TMB)) %>%
  mutate(group2 = factor(ifelse(group == "WT", "WT", "AllMutant"), levels = c("WT", "AllMutant")))

log_progress("xCell immune infiltration cohort: ",
             paste(names(table(model_df$group2)), table(model_df$group2), sep = "=", collapse = ", "))

cell_types <- setdiff(colnames(model_df), c("patient_id", "group", "TMB", "msi_score", "group2"))

run_model <- function(cell_type) {
  f <- as.formula(paste0("`", cell_type, "` ~ group2 + msi_score + TMB"))
  fit <- lm(f, data = model_df)
  cf <- summary(fit)$coefficients
  data.frame(cell_type = cell_type, beta = cf["group2AllMutant", "Estimate"],
             se = cf["group2AllMutant", "Std. Error"], p_value = cf["group2AllMutant", "Pr(>|t|)"])
}

results <- bind_rows(lapply(cell_types, run_model))
results$p_adj_BH <- p.adjust(results$p_value, method = "BH")
results <- results %>% arrange(p_value)
save_table(results, "immune_infiltration_allmutant_vs_wt_xcell", ANALYSIS)
log_progress("xCell scores, AllMutant vs WT, MSI+TMB-adjusted (", nrow(results), " cell types):\n",
             paste(capture.output(print(results, digits = 3)), collapse = "\n"))

## Cross-method comparison table: does xCell agree with quanTIseq on
## direction/significance, cell-type by cell-type where names overlap?
## Cell-type naming differs slightly between tools (e.g. xCell's finer
## subtypes), so this is a best-effort name-matched comparison, not a
## guaranteed 1:1 mapping.
qts <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, "immune_infiltration_allmutant_vs_wt_quantiseq.csv"))
comparison <- inner_join(
  results %>% transmute(cell_type, beta_xcell = beta, p_xcell = p_value, p_adj_xcell = p_adj_BH),
  qts %>% transmute(cell_type, beta_quantiseq = beta, p_quantiseq = p_value, p_adj_quantiseq = p_adj_BH),
  by = "cell_type"
)
save_table(comparison, "immune_infiltration_xcell_vs_quantiseq_comparison", ANALYSIS)
log_progress("Cell types with name-matched overlap between xCell and quanTIseq: ", nrow(comparison))
if (nrow(comparison) > 0) {
  log_progress(paste(capture.output(print(comparison, digits = 3)), collapse = "\n"))
}

## Plot: all cell types, AllMutant vs WT
plot_df <- model_df %>%
  tidyr::pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "score")
p <- ggplot(plot_df, aes(x = group2, y = score, fill = group2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.6) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 6) +
  labs(x = NULL, y = "xCell score", title = "xCell immune/stromal scores: AllMutant vs WT") +
  theme_kcna5(base_size = 8) + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(p, "immune_infiltration_allmutant_vs_wt_xcell", ANALYSIS, width = 15, height = 11)

log_progress("xCell immune infiltration (AllMutant vs WT) analysis complete.")
