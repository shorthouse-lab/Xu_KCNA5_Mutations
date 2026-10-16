# Positive selection and structural disruption of KCNA5 by mutation in colorectal cancer

Analysis code accompanying the manuscript (submitted to PLOS Computational Biology).
KCNA5 (Kv1.5), a voltage-gated potassium channel with no prior role in colorectal
cancer genomics, is recurrently mutated in TCGA-COADREAD. This repository combines
domain-resolved selection analysis (dN/dS), structural and molecular-dynamics
modelling, and transcriptomic enrichment to ask whether these mutations are
functionally meaningful, and whether their consequences depend on where in the
protein they fall.

**See `FIGURE_CODE_MANIFEST.md` for the complete figure-by-panel mapping to scripts.**
This README gives the folder-level overview; the manifest gives per-panel detail,
including a few places where the exact script used for the final figure still needs
confirming among several similar iterations.

## Data sources

- **COSMIC** (cancer.sanger.ac.uk, v103) — somatic mutations, "large intestine" /
  "carcinoma" / "adenocarcinoma" filter. Not included in this repository (licensing);
  see `dNdS/README.md` for download instructions.
- **TCGA-COADREAD** RNA-seq (raw counts) and masked somatic mutation MAFs — via the
  UCSC Xena Browser / GDC.
- **cBioPortal** — clinical annotations (MSI score, tumour mutation burden) and
  KCNA5 copy-number calls (GISTIC).
- **SWISS-MODEL** homology model of KCNA5, UniProt P22460, template 2r9r.3.B
  (residues 119–527).
- **MSigDB** (via `msigdbr`) — Hallmark, KEGG, Reactome, and Oncogenic Signatures
  (C6) gene sets.

## Folder structure

### Analyses
| Folder | Contents |
|---|---|
| `dNdS/` | Genome-wide dndscv selection analysis across potassium channel genes (Fig 1A) |
| `Clustering/` | NMC spatial clustering, whole-protein lollipop plot, selectivity-filter proximity permutation test (Fig 1D, 2A, S1) |
| `binding_ddg/` | FoldX inter-subunit binding ΔΔG, T1 domain interface, both chain orientations |
| `md/` | Raw GROMACS molecular dynamics trajectories (WT and D166N, 100 ns × 3 replicates each) — CHARMM36m force field, POPC bilayer, TIP3P water (see `md/KCNA5_WT_100ns/gromacs/`) |
| `Experimental_correlation/` | Comparison against published electrophysiology (KCNA2 alanine-scanning) |
| `patient_groups/` | Defines the T1 / S1-S4 / WT (and domain-split) patient groupings from the TCGA MAF + CNA calls, used by every transcriptomic script below | — (upstream prerequisite) |
| `genomics_selection/` | Domain-resolved dN/dS (T1 vs rest, 4-domain split, whole-gene), COSMIC mutation extraction, selectivity-filter-avoidance permutation tests | 1A (if re-derived), S1 |
| `pore_profiling/` | HOLE pore-radius and CHAP hydrophobicity profiling of the WT structure and point mutants | 2B |
| `biophysical_fingerprint/` | 23-feature biophysical characterisation of every possible pore-domain substitution — FoldX energetic decomposition, pore geometry (HOLE/CHAP), elastic-network dynamics, key H-bond distances — embedded by PCA (`44`) and ranked by per-feature effect size (`49`), compared against a literature reference panel and codon-reachable background by permutation test. Full FoldX/HOLE/CHAP batch-generation chain included | 3 |
| `t1_interface_energetics/` | FoldX binding/folding ΔΔG landscape for T1-domain interface mutations, compared against the KCNA2 alanine-scan experimental dataset | 5D, 5E |
| `gnm_elastic_network/` | Gaussian Network Model (ANM) of the WT tetramer; positions recurrently-mutated T1 residues relative to predicted dynamic hinge sites | S2 |
| `md_analysis/` | Ramachandran (backbone dihedral) analysis of the G128R MD trajectory vs WT | 6D |
| `lollipop_plots/` | Domain-zoomed lollipop plots (T1, S1-S4, pore/S6) | 4, 5A |
| `transcriptomics/genomewide_allmutant/` | DESeq2 (AllMutant vs WT, KCNA5-expression groups), the corrected GSEA pipeline, volcano plot, Hallmark pathway barplot, full pathway-universe UMAP landscape | 7A, 7B, 7C, S3 |
| `transcriptomics/domain_split_t1_s1s4/` | DESeq2 and GSEA restricted to T1-mutant and S1-S4-mutant patients; direct T1-vs-S1-S4 contrast; gene-overlap Venn; pathway concordance scatter; hypoxia validation score | 8A, 8B, 8C, 8D |
| `transcriptomics/expression_vs_mutation/` | Does KCNA5 mutation phenocopy low KCNA5 expression? Gene- and pathway-level concordance between the two | S5, S6 |
| `transcriptomics/immune_infiltration/` | xCell immune cell deconvolution, WT vs mutant tumours | S4 |
| `shared_R_helpers/` | `constants.R`, `progress.R`, `plot_theme.R` — sourced by every R script above. Point each script's `source("R/...")` calls here (or copy into your own working directory) before running standalone. |

## A note on the GSEA ranking statistic

Every GSEA result in this repository (all of Figures 7 and 8, and the Hallmark
panel of Figure S6) is computed by **`transcriptomics/genomewide_allmutant/
run_gsea_master_statranked_fix.R`**, which ranks genes by DESeq2's own Wald
`stat` column. This is deliberate: an earlier version of this pipeline ranked
genes by `sign(log2FoldChange) × -log10(padj)`, which collapses the majority of
genes into tied values given this cohort's padj distribution (only ~5,900 unique
padj values across ~28,000 tested genes). fgsea's own documentation warns that
heavily tied preranked statistics make results order-dependent. The Wald `stat` column is continuous (zero ties) and
gives a stable, reproducible result.

## Reproduction order

1. `patient_groups/define_groups.R`, `define_domain_groups.R`
2. `transcriptomics/genomewide_allmutant/run_deseq2_contrast_all_mutant.R`,
   `run_deseq2_kcna5_expression_groups.R`;
   `transcriptomics/domain_split_t1_s1s4/run_deseq2_contrasts_domain.R`
3. `transcriptomics/genomewide_allmutant/run_gsea_master_statranked_fix.R`
   (the single source of all GSEA results used downstream)
4. Any plotting/concordance script in `transcriptomics/*` — each reads the
   tables produced by steps 1–3 and will raise an explicit error naming the
   missing file if run out of order.
5. Structural analyses (`genomics_selection/`, `pore_profiling/`,
   `biophysical_fingerprint/`, `t1_interface_energetics/`,
   `gnm_elastic_network/`, `lollipop_plots/`, `md_analysis/`) are independent
   of the transcriptomics pipeline and of each other's order, but several
   depend on COSMIC mutation extraction (`genomics_selection/
   01_prepare_cosmic_mutations.R`) and the SWISS-MODEL structure/FoldX runs
   having already been generated — see each script's own header.

## Dependencies

### R
DESeq2, fgsea (≥1.26), dndscv, bio3d, GSVA, dplyr, ggplot2, ggrepel, uwot,
SummarizedExperiment, BiocParallel, msigdbr

### Python
pandas, numpy, matplotlib, seaborn, scipy, biopython, prody

### External tools
FoldX 5.1, HOLE, CHAP (Channel Annotation Package, via Docker), GROMACS 2021.5,
MODELLER, Clustal Omega, PlotDigitizer

## Citation

If you use this code, please cite Xu et al.

## Note

This folder and documentation (as well as some aspects of the scripts) was made with the help of Claude Code. 
