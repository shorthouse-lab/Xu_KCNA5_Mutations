# Xu_KCNA5_Mutations

Computational analysis of cancer-associated mutations in KCNA5 (Kv1.5) in colorectal
adenocarcinoma (COAD).

## Repository structure

```
Xu_KCNA5_Mutations/
├── dNdS/                 # dN/dS positive selection analysis (R, Python)
├── clustering/             # Spatial mutation clustering and lollipop plot (R, Python)
├── binding_ddg/            # Binding ΔΔG calculations using FoldX (Python)
├── md/                     # Molecular dynamics simulation and analysis (GROMACS, Python)
└── pathway_analysis/       # PROGENY, DESeq2, and GSEA on TCGA-COAD (R)
```

Each directory contains its own README with full instructions for data, dependencies,
and reproduction.

## Quick start

Clone the repository:

```bash
git clone https://github.com/shorthouse-lab/Xu_KCNA5_Mutations.git
cd Xu_KCNA5_Mutations
```

Then follow the README in each subdirectory. The analyses are independent of each
other and can be run in any order, with the exception of `pathway_analysis/` where
scripts must be run sequentially (see `pathway_analysis/README.md`).

## Analysis overview

### 1. dN/dS (`dNdS/`)

Identifies potassium channel genes under positive selection in COAD using COSMIC
mutation data and `dndscv`. Also includes a Monte Carlo simulation testing whether
KCNA5 mutations cluster near the selectivity filter beyond random expectation.

Tools: R (`dndscv`), Python (`scipy`, `matplotlib`)  
Data: COSMIC Genome Screens Mutant (large intestine, adenocarcinoma)

### 2. Mutation clustering (`clustering/`)

Identifies statistically significant spatial clusters of KCNA5 missense mutations
along the protein sequence using `iPAC`, and visualises them as an annotated
lollipop plot with domain boundaries and cluster density.

Tools: R (`iPAC`), Python (`matplotlib`)  
Data: COSMIC (shared with `dNdS/`)

### 3. Binding ΔΔG (`binding_ddg/`)

Calculates the effect of all possible missense mutations at the KCNA5 T1 domain
inter-subunit interface on binding energy using FoldX, across both orientations
of the interface. Results are cross-referenced against COSMIC COAD mutations.

Scripts adapted from: https://github.com/shorthouse-lab/binding_ddg  
Tools: Python (`biopython`), FoldX 5.1  
Data: SWISS-MODEL structure (UniProt P22460, residues 120–527, template 9eef.1.D)

### 4. Molecular dynamics (`md/`)

100 ns GROMACS simulations of the KCNA5 T1 domain tetramer (WT and D166N mutant,
3 replicates each). Analysis focuses on inter-subunit distances at the buried polar
interface formed by D166 and its hydrogen-bonding partners.

Scripts adapted from: https://github.com/shorthouse-lab/Xu_KCNA5_Mutations/tree/main/md  
Tools: GROMACS, Python (`MDAnalysis`, `simple_mda`)

### 5. Pathway analysis (`pathway_analysis/`)

R pipeline analysing KCNA5 expression and T1 domain mutation status in 514 TCGA-COAD
samples. Includes PROGENY pathway scoring, DESeq2 differential expression (tumour vs.
normal; T1 mutants vs. WT), and GSEA using Hallmark gene sets.

Tools: R (`DESeq2`, `progeny`, `fgsea`)  
Data: TCGA-COAD RNA-seq (GDC Hub via cBioPortal)

## Dependencies

See each subdirectory README for specific package requirements. A summary:

| Analysis        | Language | Key tools                          |
|-----------------|----------|------------------------------------|
| dN/dS           | R, Python | dndscv, scipy                     |
| Clustering      | R, Python | iPAC, matplotlib                   |
| Binding ΔΔG     | Python   | biopython, FoldX 5.1               |
| MD simulation   | Python   | GROMACS, MDAnalysis                |
| Pathway analysis| R        | DESeq2, progeny, fgsea             |

## Data availability

Raw data files are not included in this repository:
- **COSMIC**: download from https://cancer.sanger.ac.uk/cosmic/download (requires
  free registration); filter for large intestine adenocarcinoma
- **TCGA-COAD**: download from https://www.cbioportal.org (Colorectal Cancer Pan
  Atlas, TCGA); see `pathway_analysis/README.md` for exact steps
- **MD trajectories**: not included due to file size; GROMACS input files are
  provided to reproduce simulations from the equilibrated structure

## Citation
