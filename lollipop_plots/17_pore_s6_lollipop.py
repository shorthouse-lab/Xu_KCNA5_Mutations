#!/usr/bin/env python3
"""
Python lollipop-style figure for the S5 + pore/selectivity-filter +
S6/C-term region, in the same style as 16_s4_gating_charge_precedent.py's
lollipop panel: COSMIC mutation counts (reusing prepare_gene_mutations
from Xu_KCNA5_Mutations-main/Clustering/lollipop_plot.ipynb, same source
data) with the S5/pore/S6 domain boxes below the backbone.

The pore-radius and hydrophobicity precedent panels that used to sit
below this were split out into their own figure
(21_pore_radius_hydrophobicity_mutants.py) -- they're plotted against
axial position along the pore, not residue number (residue number isn't
monotonic with axial position; see that script's docstring), which
doesn't share an x-axis with this lollipop at all, and they now overlay
WT vs. mutant structures rather than being a single fixed line.

No axvspan reference-residue highlighting and no hotspot-vs-other marker
recolouring here (unlike the S4 figure) -- all mutations are plotted in a
single uniform style, per feedback that highlighting didn't add anything.
"""
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

RESULTS_DIR = "results"
ANALYSIS = "structural"
XU_REPO = "Xu_KCNA5_Mutations-main/Clustering"

S5_START, S5_END = 433, 456
PORE_START, PORE_END = 468, 510
S6_START, S6_END = 511, 613
REGION_START, REGION_END = S5_START - 6, 540


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
# B. COSMIC mutations, pore + S6/C-term span
# =============================================================================
df = prepare_gene_mutations(f"{XU_REPO}/data/cosmic_LI_AC_filtered.tsv", "KCNA5")
region_df = df[(df["Residue"] >= REGION_START) & (df["Residue"] <= REGION_END)].copy()

DOMAIN_COLOURS = {"S5": "#FF8D00", "Pore/selectivity filter": "#808080", "S6/C-term": "#FF8D00"}

# =============================================================================
# C. Figure: lollipop only
# =============================================================================
fig, ax_lollipop = plt.subplots(figsize=(12, 5))

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

for name, start, end in [("S5", S5_START, S5_END), ("Pore/selectivity filter", PORE_START, PORE_END), ("S6/C-term", S6_START, S6_END)]:
    rect = patches.Rectangle((start, box_bottom), end - start, box_height,
                              facecolor=DOMAIN_COLOURS[name], edgecolor="k", linewidth=2, zorder=6, clip_on=True)
    ax_lollipop.add_patch(rect)
    visible_start, visible_end = max(start, REGION_START), min(end, REGION_END)
    ax_lollipop.text((visible_start + visible_end) / 2, backbone_y, name, ha="center", va="center",
                      fontsize=11, weight="bold", color="white", zorder=7, clip_on=True)

for xi, yi, lab in zip(x, y, labels):
    ax_lollipop.text(xi, yi + 0.1, lab, ha="left", va="bottom", fontsize=10, rotation=45,
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

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/pore_s6_mutations_python.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/pore_s6_mutations_python.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
