#!/usr/bin/env python3
"""
Zoomed lollipop for the T1 domain (120-216), same style as
17_pore_s6_lollipop.py's S5/pore/S6 figure -- single domain box, uniform
marker colour, no additional highlighting -- just a different region of the
sequence.
"""
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

RESULTS_DIR = "results"
ANALYSIS = "structural"
XU_REPO = "Xu_KCNA5_Mutations-main/Clustering"

## T1_START moved from 128 to 120 (this figure only -- the project-wide
## T1_DOMAIN_START/END constants in R/constants.R stay at 128-216, since the
## TCGA cohort that drives every DESeq2/GSEA/genomics domain-split analysis
## has zero mutations in 120-127 and is therefore unaffected either way).
## Literature basis: UniProt P22460's own Tetramerization-domain annotation
## (REGION 1-211, citing Burg et al. 2010) is far wider than either 120 or
## 128; this project's own domain table elsewhere (e.g. 19/20_*_pore_
## radius.py) separately calls 1-108/1-120 "Disordered (N-term)" -- so 120
## is picked as the boundary where that disordered region ends, folding the
## former 121-127 "Linker" carve-out into T1 itself. This also moves R121
## (COSMIC R121H/R121C) from just-outside the domain box to inside it.
T1_START, T1_END = 120, 216
REGION_START, REGION_END = 110, T1_END + 6


# =============================================================================
# A. Reused verbatim from lollipop_plot.ipynb
# =============================================================================
def prepare_gene_mutations(filepath, gene_symbol, sep="\t"):
    """Load mutation data, filter missense variants for a given gene,
    parse amino acid changes, and summarize counts per residue."""
    df = pd.read_csv(filepath, sep=sep)
    df = df[(df["GENE_SYMBOL"] == gene_symbol) & (df["MUTATION_DESCRIPTION"] == "missense_variant")].copy()

    residues, mutations = [], []
    for aa_change in df["MUTATION_AA"]:
        if not isinstance(aa_change, str) or not aa_change.startswith("p."):
            continue
        mut = aa_change[2:]
        if mut in {"?", ""} or "*" in mut:
            continue
        if mut[0] == mut[-1]:
            continue
        residue_match = re.search(r"\d+", mut)
        if not residue_match:
            continue
        residues.append(int(residue_match.group()))
        mutations.append(mut)

    if not residues:
        return pd.DataFrame(columns=["Residue", "Mutations", "Count"])

    data = pd.DataFrame({"Residue": residues, "Mutation": mutations})
    summary = (
        data.groupby("Residue")
        .agg(Mutations=("Mutation", lambda x: "/".join(sorted(set(x)))), Count=("Mutation", "count"))
        .reset_index()
    )
    return summary


# =============================================================================
# B. COSMIC mutations, T1 domain span
# =============================================================================
df = prepare_gene_mutations(f"{XU_REPO}/data/cosmic_LI_AC_filtered.tsv", "KCNA5")
region_df = df[(df["Residue"] >= REGION_START) & (df["Residue"] <= REGION_END)].copy()

DOMAIN_COLOURS = {"T1": "#FF8D00"}

# =============================================================================
# C. Figure: lollipop only
# =============================================================================
fig, ax_lollipop = plt.subplots(figsize=(10, 5))

x = region_df["Residue"].tolist()
y = region_df["Count"].tolist()
labels = region_df["Mutations"].tolist()

backbone_y = 0.0
box_height = 0.6
box_bottom = backbone_y - box_height / 2
ax_lollipop.hlines(y=backbone_y, xmin=REGION_START, xmax=REGION_END, linewidth=6, color="#808080", zorder=1)
for xi, yi in zip(x, y):
    ax_lollipop.plot([xi, xi], [0, yi], color="k", linewidth=1.5, zorder=3)
ax_lollipop.scatter(x, y, s=110, c="#1f77b4", edgecolors="k", linewidths=1, zorder=4)

for name, start, end in [("T1", T1_START, T1_END)]:
    rect = patches.Rectangle((start, box_bottom), end - start, box_height,
                              facecolor=DOMAIN_COLOURS[name], edgecolor="k", linewidth=2, zorder=6, clip_on=True)
    ax_lollipop.add_patch(rect)
    visible_start, visible_end = max(start, REGION_START), min(end, REGION_END)
    ax_lollipop.text((visible_start + visible_end) / 2, backbone_y, name, ha="center", va="center",
                      fontsize=11, weight="bold", color="white", zorder=7, clip_on=True)

## V195M/F196L sit one residue apart, close enough that their labels
## collide at this width -- nudge just these two slightly apart in x.
LABEL_NUDGE = {"V195M": -0.6, "F196L": 0.6}

for xi, yi, lab in zip(x, y, labels):
    label_x = xi + LABEL_NUDGE.get(lab, 0)
    ax_lollipop.text(label_x, yi + 0.1, lab, ha="left", va="bottom", fontsize=10, rotation=45,
                      color="black", zorder=5)

ax_lollipop.set_xlim(REGION_START, REGION_END)
ax_lollipop.set_ylim(box_bottom - 0.5, max(y) + 1.8)
ax_lollipop.set_yticks(range(0, max(y) + 1))
ax_lollipop.set_ylabel("Number of\nMutations", fontsize=14)
ax_lollipop.set_xlabel("KCNA5 residue", fontsize=14)
ax_lollipop.tick_params(axis="both", width=2, labelsize=13)
for spine in ax_lollipop.spines.values():
    spine.set_linewidth(2)
sns.despine(ax=ax_lollipop)

fig.tight_layout()

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/t1_interface_lollipop.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/t1_interface_lollipop.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
