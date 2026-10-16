## ---------------------------------------------------------------------------
## Gene-level companions to the pathway-level domain comparison
## (run_domain_signature_comparison.R, run_pathway_comparison_dumbbell.R):
##   A. Volcano plots, T1_vs_WT and S1S4_vs_WT, same axis limits/thresholds
##      so the two are visually comparable at a glance.
##   B. Venn diagram of significant genes (padj<0.05) shared between the two
##      contrasts -- the gene-level counterpart to the pathway dumbbell's
##      "14 T1-specific, 1 shared, 0 S1S4-specific" finding.
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(eulerr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_domain_split"
PADJ_THRESHOLD <- 0.05

t1   <- read.csv(file.path(RESULTS_DIR, "tables", "analysis3_transcriptomic", "deseq2_T1_vs_WT.csv")) %>%
  mutate(contrast = "T1 vs WT")
s1s4 <- read.csv(file.path(RESULTS_DIR, "tables", "analysis3_domain_split", "deseq2_S1S4_vs_WT_domain.csv")) %>%
  mutate(contrast = "S1S4 vs WT")

## =============================================================================
## A. Volcano plots -- shared axis limits for direct visual comparison
## =============================================================================
both <- bind_rows(t1, s1s4) %>% filter(!is.na(padj))
xlim_shared <- max(abs(both$log2FoldChange), na.rm = TRUE) * c(-1, 1)
ylim_shared <- c(0, max(-log10(both$pvalue), na.rm = TRUE) * 1.05)

make_volcano <- function(df, title_text, n_label = 15) {
  df <- df %>%
    mutate(sig = padj < PADJ_THRESHOLD,
           neg_log10_p = -log10(pvalue))
  n_sig_up <- sum(df$sig & df$log2FoldChange > 0, na.rm = TRUE)
  n_sig_down <- sum(df$sig & df$log2FoldChange < 0, na.rm = TRUE)

  label_df <- df %>% filter(sig) %>% slice_min(padj, n = n_label)

  ggplot(df, aes(x = log2FoldChange, y = neg_log10_p)) +
    geom_point(data = df %>% filter(!sig), colour = "grey75", size = 0.8, alpha = 0.6) +
    geom_point(data = df %>% filter(sig), aes(colour = log2FoldChange > 0), size = 1.4, alpha = 0.8) +
    geom_text_repel(data = label_df, aes(label = gene_symbol), size = 3, max.overlaps = 25,
                     segment.size = 0.3, min.segment.length = 0) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
    scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"), guide = "none") +
    coord_cartesian(xlim = xlim_shared, ylim = ylim_shared) +
    labs(x = "log2 fold change", y = expression(-log[10](p)),
         title = title_text,
         subtitle = paste0(n_sig_up, " up / ", n_sig_down, " down at padj<", PADJ_THRESHOLD,
                            " (n=", nrow(df), " genes tested)")) +
    theme_kcna5(base_size = 16)
}

p_volcano_t1   <- make_volcano(t1, "T1-mutant vs WT")
p_volcano_s1s4 <- make_volcano(s1s4, "S1S4-mutant vs WT")

save_fig(p_volcano_t1, "volcano_T1_vs_WT", ANALYSIS, width = 8, height = 7)
save_fig(p_volcano_s1s4, "volcano_S1S4_vs_WT", ANALYSIS, width = 8, height = 7)

log_progress("Volcano plots saved (shared axis limits x=", round(xlim_shared[2], 2),
             ", y=", round(ylim_shared[2], 2), ")")

## =============================================================================
## B. Venn diagram of significant genes (padj<0.05), T1 vs WT vs S1S4 vs WT
## =============================================================================
t1_sig   <- t1 %>% filter(!is.na(padj), padj < PADJ_THRESHOLD) %>% pull(ensembl_id)
s1s4_sig <- s1s4 %>% filter(!is.na(padj), padj < PADJ_THRESHOLD) %>% pull(ensembl_id)

log_progress("Significant genes: T1_vs_WT=", length(t1_sig), ", S1S4_vs_WT=", length(s1s4_sig),
             ", shared=", length(intersect(t1_sig, s1s4_sig)))

venn_fit <- euler(list(`T1 vs WT` = t1_sig, `S1S4 vs WT` = s1s4_sig))

png(file.path(RESULTS_DIR, "figures", ANALYSIS, "venn_sig_genes_T1_vs_S1S4.png"),
    width = 7, height = 6, units = "in", res = FIG_DPI)
print(plot(venn_fit, quantities = TRUE, fills = c("#d95f02", "#7570b3"), alpha = 0.6,
           labels = list(fontsize = 14), quantities_fontsize = 16,
           main = paste0("Significant genes (padj<", PADJ_THRESHOLD, "): T1 vs WT vs S1S4 vs WT")))
dev.off()

pdf(file.path(RESULTS_DIR, "figures", ANALYSIS, "venn_sig_genes_T1_vs_S1S4.pdf"), width = 7, height = 6)
print(plot(venn_fit, quantities = TRUE, fills = c("#d95f02", "#7570b3"), alpha = 0.6,
           labels = list(fontsize = 14), quantities_fontsize = 16,
           main = paste0("Significant genes (padj<", PADJ_THRESHOLD, "): T1 vs WT vs S1S4 vs WT")))
dev.off()

venn_counts <- data.frame(
  category = c("T1-specific", "S1S4-specific", "Shared"),
  n_genes = c(length(setdiff(t1_sig, s1s4_sig)), length(setdiff(s1s4_sig, t1_sig)), length(intersect(t1_sig, s1s4_sig)))
)
save_table(venn_counts, "venn_sig_genes_T1_vs_S1S4_counts", ANALYSIS)
log_progress("Venn diagram complete:\n", paste(capture.output(print(venn_counts)), collapse = "\n"))
