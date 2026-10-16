# Figure-to-code manifest

Maps every figure/panel in `Xuetal_KCNA5_colorectalcancer_v3.docx` to the script(s)
that produced it. Built by tracing each figure caption against the analysis
pipeline in `Analysis/scripts/`. Scope: **only code for figures actually in the
manuscript** (main Figures 1-8 + Figures S1-S6), not the full exploratory pipeline.

Status key: ✅ copied in, confirmed · ⚠️ copied in, best guess among several
similar candidate scripts, not yet narrowed to one · ❓ not located in this
pipeline.

**Confirmed out of scope (2026-10-16):** every remaining ❓ item (Fig 1B/C/D,
2A, 6A/B/C/E) is a manual structure render or an MD-trajectory analysis
already covered by Binglin's own folders (`Clustering/`, `md/`,
`Experimental_correlation/`) — no further action needed. The ⚠️ items in
Figures 2 and 5 (multiple similarly-named script iterations) are left as-is
for you to pick the exact matching one during curation; not investigating
further.

All R scripts `source("R/constants.R")`, `source("R/progress.R")`,
`source("R/plot_theme.R")` — these are copied once into `shared_R_helpers/`;
point the `source()` calls at that path (or copy into each script's own
working directory) before running anything here standalone.

## Figure 1 — KCNA5 is under mutational selection in CRC
| Panel | Description | Script | Status |
|---|---|---|---|
| A | dN/dS of potassium channels expressed in colorectal tract (multi-gene, not KCNA5-domain-split) | `dNdS/cosmic_dnds_calculation.R` (Binglin's original) | ❓ confirm this is still the source — Methods text matches Binglin's original multi-gene dndscv run, not any domain-split script I built |
| B, C | Structure renders (full tetramer; single subunit, domains highlighted) | — | ❓ likely a PyMOL/ChimeraX session, not a script in this pipeline |
| D | Whole-protein lollipop plot of observed mutations | `Clustering/lollipop_plot.ipynb` (Binglin's original) | ❓ confirm — no equivalent found in the newer pipeline |

## Figure 2 — Mutations avoid the selectivity filter
| Panel | Description | Script | Status |
|---|---|---|---|
| A | Lollipop with selectivity-filter-avoidance region highlighted | — | ❓ not located as a standalone script; may be a variant of Binglin's lollipop notebook |
| B | HOLE pore-profile render, mutations highlighted | `pore_profiling/20_hole_pore_radius.py` or `19_direct_pore_radius.py` | ⚠️ both copied, pick the one matching the actual rendered figure |
| S1 | Permutation test for selectivity-filter mutation avoidance (two Monte Carlo schemes) | `genomics_selection/04b_permutation_nn_distance.R`, `04c_permutation_nn_distance_tetramer.R` | ⚠️ both copied; confirm which (or both) match the two reported p-values (p=0.0213 full-sequence, p=0.0001 pore-restricted) |

## Figure 3 — Biophysical fingerprint of pore/filter mutations
| Panel | Description | Script | Status |
|---|---|---|---|
| A | PCA of 23-feature biophysical fingerprint, pore/filter region (433-527), vs literature reference classes and codon background | `biophysical_fingerprint/44_full_landscape_pca_umap.R` | ✅ **confirmed — manuscript text now states 23 features, matching this script (not the earlier 21-feature `36`/`31` iterations, which have been removed from this folder as superseded)** |
| B | Per-feature effect size (standardised mean difference vs background), 4 mutational classes, Welch's t-test significance stars | `biophysical_fingerprint/49_feature_importance.R` | ✅ (reads `44`'s output directly) |

Full prerequisite chain for `44`/`49`, all now copied into this folder:
`24_mutation_effect_summary.py` (observed-13 mutation summary) → `32_observed_foldx_breakdown.py`
(FoldX decomposition, observed mutations) · `29_pore_background_enumeration.R` /
`37_pore_full_saturation_enumeration.R` (codon-reachable / full-saturation substitution
enumeration, STRUCT_FEATURES) · `30_pore_background_foldx.py` / `38_pore_full_saturation_foldx.py`
(FoldX ΔΔG, background) · `33_pore_background_hole.py` / `41_pore_full_saturation_hole.py` (HOLE
pore radius, background) · `34_pore_background_chap.py` / `42_pore_full_saturation_chap.py` (CHAP
hydrophobicity, background) · `39_regenerate_pdbs_for_hbonds.py` → `40_hbond_distances.py`
(H-bond distance features) · `43_reference_mutation_extraction.R` (35-mutation/5-category
literature reference set) · `23_enm_scaffold_filter_coupling.R` (ANM flexibility/filter-coupling
features, cached and reused by the enumeration scripts above).

## Figure 4 — S1-S4 domain lollipop
| Panel | Description | Script | Status |
|---|---|---|---|
| (single) | Lollipop of S1-S4 mutations | `lollipop_plots/67_s1s4_lollipop_simple.py` | ✅ |

## Figure 5 — T1 domain interface binding/folding
| Panel | Description | Script | Status |
|---|---|---|---|
| A | T1 cohort lollipop | `lollipop_plots/57_t1_interface_lollipop.py` | ✅ (note: boundary now 120-216, not 128-216 — see script header) |
| B, C | T1 structure / interface contact renders | `t1_interface_energetics/02_extract_interface_residues.R` | ⚠️ data source for contacts; the render itself is likely manual (PyMOL) |
| D | Binding/folding ΔΔG plane, cohort vs KCNA2 alanine-scan experimental comparison | `t1_interface_energetics/61_experimental_binding_vs_folding.py` or `63_experimental_binding_vs_folding_marginals.py` | ⚠️ several similarly-named iterations exist (`28`, `54`, `60`-`64`) — all copied, **pick the one matching the final figure** |
| E | Binding/folding ΔΔG plane, own cohort (G128R/R121/R167) | `t1_interface_energetics/64_observed_binding_vs_folding_samepos_marginals.py` or `60_binding_vs_folding_quadrant.py` | ⚠️ same ambiguity as D |

## Figure S2 — GNM hinge analysis
| Panel | Description | Script | Status |
|---|---|---|---|
| (single) | R121/G128/D166/R167 vs elastic-network hinge sites | `gnm_elastic_network/66_t1_gnm_hinges_mobility.py` | ✅ (boundary updated to 120-216) |

## Figure 6 — MD simulations of T1 mutations
| Panel | Description | Script | Status |
|---|---|---|---|
| A | D166 position/structure render | — | ❓ likely manual render |
| B, C | Cross-correlation of interface residue motion, WT vs D166N MD | — | ❓ **not located in this pipeline — likely Binglin's own analysis of `md/KCNA5_WT_100ns` / `md/KCNA5_D166N_100ns` trajectories, not yet in any script I have access to** |
| D | Ramachandran plot, G128 WT vs G128R | `md_analysis/71_g128r_ramachandran.py` | ✅ |
| E | G128R MD screenshot, R128/D166 proximity | — | ❓ likely a manual PyMOL/VMD screenshot, not a reproducible script |

## Figure 7 — Transcriptomics, AllMutant vs WT
| Panel | Description | Script | Status |
|---|---|---|---|
| A | Volcano, AllMutant vs WT | `transcriptomics/genomewide_allmutant/run_volcano_all_and_subclasses.R` | ✅ |
| B | Hallmark pathway barplot, AllMutant vs WT | `transcriptomics/genomewide_allmutant/24_allmutant_pathway_barplot.py` | ✅ |
| C | UMAP pathway landscape, all 4 collections, FDR q≤0.25 | `transcriptomics/genomewide_allmutant/run_pathway_umap_landscape.R` | ✅ |
| S3 | Per-collection AllMutant barplots (KEGG/Reactome/Oncogenic) | `transcriptomics/genomewide_allmutant/24_allmutant_pathway_barplot.py` (same script, other outputs) | ✅ |
| S4 | xCell immune infiltration, WT vs mutant | `transcriptomics/immune_infiltration/run_immune_infiltration_xcell.R` | ✅ |
| S5 | Gene-level expression-vs-mutation concordance scatter | `transcriptomics/expression_vs_mutation/run_kcna5_expression_vs_mutation_comparison.R` | ✅ |
| S6 | Pathway-level expression-vs-mutation concordance (56-pathway set) | `transcriptomics/expression_vs_mutation/run_kcna5_expression_vs_mutation_all_collections.R` + `run_kcna5_expression_vs_mutation_dumbbells_all.R` | ✅ |

✅ **Resolved.** All four scripts that previously recomputed GSEA internally with the old tie-heavy ranking (`run_gsea_domain.R`, `run_gsea_t1_vs_s1s4_direct_combined.R`, `run_deseq2_kcna5_expression_groups.R`, `run_kcna5_expression_vs_mutation_all_collections.R`) have been edited to remove that logic entirely — they now only *read* the `fgsea_*.csv` tables that `run_gsea_master_statranked_fix.R` produces, and will `stop()` with a clear message naming the missing file if run before it. GSEA computation lives in exactly one place in this repository now. The manuscript text has also been updated to match (confirmed by the user).

## Figure 8 — T1 vs S1-S4 domain comparison
| Panel | Description | Script | Status |
|---|---|---|---|
| A | Venn, DE genes T1 vs WT and S1-S4 vs WT | `transcriptomics/domain_split_t1_s1s4/run_volcano_venn_domain.R` | ✅ |
| B | Direct T1-vs-S1-S4 Hallmark barplot | `transcriptomics/domain_split_t1_s1s4/25_t1_vs_s1s4_direct_hallmark_barplot.py` (needs `run_gsea_t1_vs_s1s4_direct_combined.R` run first) | ✅ |
| C | Hypoxia score, T1 vs S1-S4 | `transcriptomics/domain_split_t1_s1s4/26_hypoxia_single_method_seaborn.py` | ✅ |
| D | Pathway comparison scatter, all 4 gene sets, T1 vs S1-S4 | `transcriptomics/domain_split_t1_s1s4/run_domain_signature_comparison.R` (pooled 4-collection output) | ✅ |

✅ **Resolved** — same fix as Figure 7 (see above): `run_gsea_t1_vs_s1s4_direct_combined.R` no longer computes GSEA itself, only pools the already-correct per-collection tables.

## Upstream prerequisites (not figures themselves, but required inputs)
- `patient_groups/define_groups.R`, `define_domain_groups.R` — defines the T1/S1-S4/WT patient groupings used by every transcriptomic script above.
- `genomics_selection/01_prepare_cosmic_mutations.R` — COSMIC mutation extraction, feeds the dN/dS and structural scripts.
- `transcriptomics/genomewide_allmutant/run_deseq2_contrast_all_mutant.R`, `run_deseq2_kcna5_expression_groups.R`, `transcriptomics/domain_split_t1_s1s4/run_deseq2_contrasts_domain.R` — DESeq2 differential expression, prerequisite for every volcano/GSEA script.