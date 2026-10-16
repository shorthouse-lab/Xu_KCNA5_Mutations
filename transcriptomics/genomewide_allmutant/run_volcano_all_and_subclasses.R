## ---------------------------------------------------------------------------
## Two headline volcano figures:
##   A. AllMutant_vs_WT (all 26 missense KCNA5 mutants pooled, T1+Non-T1) --
##      single panel.
##   B. The two major structural subclasses, T1_vs_WT and S1S4_vs_WT, as one
##      combined two-panel figure (shared axis limits, patchwork) rather than
##      two separate files -- the gene-level counterpart to the pathway
##      dumbbell figure.
##
## CAVEAT carried over from the Other_domain_mutant outlier check earlier in
## this session: TCGA-A6-5656 (T39M, part of the Non-T1 pool) and
## TCGA-G4-6628 (R234H, same pool) each have extreme Cook's-distance outlier
## genes (DNAH7/CCAT2 and LINC01426 respectively) that reappear in panel A's
## top-hit list diluted but not removed. These specific genes are flagged in
## the saved table (outlier_flag column) -- NOT removed from the plot, since
## that would hide real cohort composition, but should not be read as a
## genuine AllMutant transcriptional effect.
## ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "analysis3_transcriptomic"
PADJ_THRESHOLD <- 0.05
OUTLIER_GENES <- c("DNAH7", "CCAT2", "LINC01426")  # see header caveat

all_mut <- read.csv(file.path(RESULTS_DIR, "tables", "analysis3_transcriptomic", "deseq2_AllMutant_vs_WT.csv"))
t1      <- read.csv(file.path(RESULTS_DIR, "tables", "analysis3_transcriptomic", "deseq2_T1_vs_WT.csv"))
s1s4    <- read.csv(file.path(RESULTS_DIR, "tables", "analysis3_domain_split", "deseq2_S1S4_vs_WT_domain.csv"))

make_volcano <- function(df, title_text, xlim_shared, ylim_shared, n_label = 15, flag_genes = character(0)) {
  df <- df %>%
    mutate(sig = padj < PADJ_THRESHOLD, neg_log10_p = -log10(pvalue),
           is_flagged = gene_symbol %in% flag_genes)
  n_sig_up <- sum(df$sig & df$log2FoldChange > 0, na.rm = TRUE)
  n_sig_down <- sum(df$sig & df$log2FoldChange < 0, na.rm = TRUE)
  label_df <- df %>% filter(sig) %>% slice_min(padj, n = n_label)

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_p)) +
    geom_point(data = df %>% filter(!sig), colour = "grey75", size = 0.8, alpha = 0.6) +
    geom_point(data = df %>% filter(sig & !is_flagged), aes(colour = log2FoldChange > 0), size = 1.4, alpha = 0.8) +
    geom_point(data = df %>% filter(sig & is_flagged), shape = 4, colour = "black", size = 2, stroke = 1) +
    geom_text_repel(data = label_df, aes(label = gene_symbol), size = 5, max.overlaps = 25,
                     segment.size = 0.3, min.segment.length = 0) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
    scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"), guide = "none") +
    coord_cartesian(xlim = xlim_shared, ylim = ylim_shared) +
    labs(x = "log2 fold change", y = expression(-log[10](p)),
         title = title_text,
         subtitle = paste0(n_sig_up, " up / ", n_sig_down, " down at padj<", PADJ_THRESHOLD,
                            " (n=", nrow(df), " genes tested)")) +
    theme_kcna5(base_size = 16)
  p
}

## =============================================================================
## A. AllMutant vs WT -- single panel
## =============================================================================
all_mut_f <- all_mut %>% filter(!is.na(padj))
xlim_a <- max(abs(all_mut_f$log2FoldChange), na.rm = TRUE) * c(-1, 1)
ylim_a <- c(0, max(-log10(all_mut_f$pvalue), na.rm = TRUE) * 1.05)

p_all <- make_volcano(all_mut_f, "All KCNA5 missense mutants vs WT", xlim_a, ylim_a)

save_fig(p_all, "volcano_AllMutant_vs_WT", ANALYSIS, width = 8.5, height = 7.5)

## =============================================================================
## B. T1 vs WT and S1S4 vs WT -- one combined two-panel figure, shared axes
## =============================================================================
both <- bind_rows(t1 %>% filter(!is.na(padj)), s1s4 %>% filter(!is.na(padj)))
xlim_b <- max(abs(both$log2FoldChange), na.rm = TRUE) * c(-1, 1)
ylim_b <- c(0, max(-log10(both$pvalue), na.rm = TRUE) * 1.05)

p_t1   <- make_volcano(t1 %>% filter(!is.na(padj)), "T1-mutant vs WT", xlim_b, ylim_b)
p_s1s4 <- make_volcano(s1s4 %>% filter(!is.na(padj)), "S1S4-mutant vs WT", xlim_b, ylim_b)

p_combined <- p_t1 + p_s1s4 +
  plot_annotation(title = "The two major KCNA5 mutation subclasses vs WT",
                   theme = theme(plot.title = element_text(face = "bold", size = 18)))

save_fig(p_combined, "volcano_T1_and_S1S4_vs_WT_combined", ANALYSIS, width = 14, height = 7.5)

log_progress("Both volcano figures saved.")
