## ---------------------------------------------------------------------------
## Separate dumbbell figures per pathway database (Hallmark, KEGG, Reactome,
## oncogenic signatures) for the KCNA5-low-vs-high (WT) vs AllMutant-vs-WT
## concordance comparison -- one panel per collection rather than the single
## Hallmark-only dumbbell from run_kcna5_expression_vs_mutation_comparison.R,
## since the full picture (42/56 concordant, with the 14 discordant pathways
## concentrated almost entirely in cell-cycle/proliferation terms) only
## shows up once Reactome/oncogenic/KEGG are included alongside Hallmark.
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_transcriptomic"

collections <- list(
  hallmark  = list(file = "kcna5_expression_vs_mutation_hallmark_concordance.csv",  prefix = "HALLMARK_", label = "Hallmark"),
  kegg      = list(file = "kcna5_expression_vs_mutation_kegg_concordance.csv",      prefix = "KEGG_",     label = "KEGG"),
  reactome  = list(file = "kcna5_expression_vs_mutation_reactome_concordance.csv",  prefix = "REACTOME_", label = "Reactome"),
  oncogenic = list(file = "kcna5_expression_vs_mutation_oncogenic_concordance.csv", prefix = "",           label = "Oncogenic signatures")
)

for (coll_name in names(collections)) {
  spec <- collections[[coll_name]]
  pathway_df <- read.csv(file.path(RESULTS_DIR, "tables", ANALYSIS, spec$file)) %>%
    mutate(pathway_label = gsub(spec$prefix, "", pathway), pathway_label = gsub("_", " ", pathway_label),
           pathway_label = stringr::str_wrap(pathway_label, width = 40),
           concordant = sign(NES_expr) == sign(NES_mut)) %>%
    arrange(NES_mut) %>%
    mutate(pathway_label = factor(pathway_label, levels = unique(pathway_label)))

  n_path <- nrow(pathway_df)
  n_conc <- sum(pathway_df$concordant)
  log_progress(coll_name, ": ", n_path, " pathways, ", n_conc, " concordant, ", n_path - n_conc, " discordant")

  long <- bind_rows(
    pathway_df %>% transmute(pathway_label, NES = NES_expr, padj = padj_expr, contrast = "KCNA5-low vs -high (WT)"),
    pathway_df %>% transmute(pathway_label, NES = NES_mut, padj = padj_mut, contrast = "AllMutant vs WT")
  ) %>%
    mutate(contrast = factor(contrast, levels = c("KCNA5-low vs -high (WT)", "AllMutant vs WT")))

  fig_height <- max(3, 0.4 * n_path + 2)

  p <- ggplot() +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_segment(data = pathway_df, aes(x = NES_expr, xend = NES_mut, y = pathway_label, yend = pathway_label,
                                         colour = concordant), linewidth = 0.8) +
    geom_point(data = long, aes(x = NES, y = pathway_label, shape = contrast), size = 3, colour = "black") +
    scale_colour_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02"), guide = "none") +
    scale_shape_manual(values = c("KCNA5-low vs -high (WT)" = 1, "AllMutant vs WT" = 16), name = NULL) +
    labs(x = "Normalised enrichment score (NES)", y = NULL,
         title = paste0(spec$label, ": mutation vs low expression, ", n_conc, "/", n_path, " concordant"),
         subtitle = "Orange connector = opposite direction between contrasts; green = same direction") +
    theme_kcna5(base_size = 15) +
    theme(legend.position = "top", axis.text.y = element_text(size = 12, lineheight = 0.85),
          plot.title = element_text(size = 16))

  out_name <- paste0("kcna5_expression_vs_mutation_", coll_name, "_dumbbell")
  save_fig(p, out_name, ANALYSIS, width = 13, height = fig_height)
}

log_progress("Per-collection dumbbell figures complete.")
