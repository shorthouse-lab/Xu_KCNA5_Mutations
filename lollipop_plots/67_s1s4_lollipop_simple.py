#!/usr/bin/env python3
"""
Zoomed lollipop for the full S1-S4 voltage-sensor domain (247-416), same
simple style as 17_pore_s6_lollipop.py / 57_t1_interface_lollipop.py /
66_s4_lollipop_simple.py -- uniform marker colour, no extra highlighting --
just extended to cover all four VSD segments instead of S4 alone, matching
the domain boundaries used in 15_vsd_charged_residue_enrichment.R
(S1 247-270, S2 322-344, S3 345-376, S4 392-416; the S1-S2 and S3-S4 spans
are linker, not drawn as a domain box, same convention as the pore/S6
lollipop's gap between S5/pore/S6).
"""
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

RESULTS_DIR = "results"
ANALYSIS = "structural"
XU_REPO = "Xu_KCNA5_Mutations-main/Clustering"

SEGMENTS = [("S1", 247, 270), ("S2", 322, 344), ("S3", 345, 376), ("S4", 392, 416)]
REGION_START, REGION_END = SEGMENTS[0][1] - 6, SEGMENTS[-1][2] + 6


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
# B. COSMIC mutations, S1-S4 span
# =============================================================================
df = prepare_gene_mutations(f"{XU_REPO}/data/cosmic_LI_AC_filtered.tsv", "KCNA5")
region_df = df[(df["Residue"] >= REGION_START) & (df["Residue"] <= REGION_END)].copy()

## Matches the master domain palette (Xu_KCNA5_Mutations-main/Clustering/
## lollipop_plot.ipynb's KCNA5_domains list, and 16_s4_gating_charge_
## precedent.py's DOMAIN_COLOUR for S4) -- S1-S4 are all the same red there,
## not individually distinguished.
DOMAIN_COLOURS = {"S1": "#FF0D0D", "S2": "#FF0D0D", "S3": "#FF0D0D", "S4": "#FF0D0D"}

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

for name, start, end in SEGMENTS:
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

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/s1s4_lollipop_simple.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/s1s4_lollipop_simple.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
