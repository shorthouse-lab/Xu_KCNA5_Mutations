#!/usr/bin/env python3
"""
Merges 61_experimental_binding_vs_folding.py (the 30-mutation KCNA2
alanine-scan reference set, coloured by experimental deltaV1/2) with
60_binding_vs_folding_quadrant.py's marginal-histogram layout (each axis's
basal/background distribution shown directly alongside the scatter, rather
than an arbitrary threshold line) -- so the experimental-colour scatter can
also be read against the full shape of what's biophysically possible at
these positions.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from adjustText import adjust_text

RESULTS_DIR = "results"
ANALYSIS = "structural"
BG_COLOR = "#2a78d6"

exp = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/experimental_binding_vs_folding.csv")
folding_bg = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_domain_folding_saturation_ddg_progress.csv")["ddG_kcal_mol"]
binding_bg_df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/binding_ddg_expected_vs_observed.csv")
binding_bg = binding_bg_df.loc[binding_bg_df["Group"] == "Expected", "DDG"]

usable = exp[~exp["unmeasurable"] & exp["max_abs_binding_ddg"].notna() & exp["folding_ddG"].notna()].copy()
unmeasurable = exp[exp["unmeasurable"] & exp["max_abs_binding_ddg"].notna() & exp["folding_ddG"].notna()].copy()

fig = plt.figure(figsize=(13, 11))
gs = GridSpec(2, 3, width_ratios=(5, 1.2, 0.3), height_ratios=(1.2, 5),
              wspace=0.08, hspace=0.05, figure=fig)

ax_scatter = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatter)
cax = fig.add_subplot(gs[1, 2])

ax_scatter.scatter(unmeasurable["max_abs_binding_ddg"], unmeasurable["folding_ddG"], facecolors="none",
                    edgecolors="#555555", linewidths=1.5, marker="s", s=100, zorder=3,
                    label="Expression-null (unmeasurable)")
sc = ax_scatter.scatter(usable["max_abs_binding_ddg"], usable["folding_ddG"], c=usable["x"],
                         cmap="coolwarm", vmin=-15, vmax=15, s=110, edgecolors="k", linewidths=0.8, zorder=4)

x_all = pd.concat([usable["max_abs_binding_ddg"], unmeasurable["max_abs_binding_ddg"]])
y_all = pd.concat([usable["folding_ddG"], unmeasurable["folding_ddG"]])

## Anchor each label exactly at its point's true (x, y) -- va="bottom"
## alone already renders the text just above that anchor, so the marker
## doesn't need a manually-offset y (adjust_text draws its leader line back
## to this anchor, so offsetting it here made the line visibly miss the
## actual point).
texts = [ax_scatter.text(row.max_abs_binding_ddg, row.folding_ddG, row.kcna5_label,
                          fontsize=13, fontweight="bold", ha="center", va="bottom")
         for row in pd.concat([usable, unmeasurable]).itertuples()]
## Kept deliberately mild -- a stronger expand/force_text spreads labels
## further apart, but past a point that means shoving some labels well
## away from their point just to avoid a neighbour, producing long
## crossing leader lines. Mild settings tolerate a little closeness/
## overlap in the most crowded cluster in exchange for every line staying
## short.
adjust_text(texts, ax=ax_scatter, arrowprops=dict(arrowstyle="-", color="grey", lw=0.6),
            expand=(1.15, 1.25), force_text=(0.25, 0.3))

ax_scatter.axvline(0, color="#999999", lw=1, zorder=1)
ax_scatter.axhline(0, color="#999999", lw=1, zorder=1)
ax_scatter.set_xlabel("Binding ΔΔG (kcal/mol)", fontsize=12)
ax_scatter.set_ylabel("Folding ΔΔG (kcal/mol)", fontsize=12)
ax_scatter.tick_params(axis="both", labelsize=13)
ax_scatter.legend(loc="lower right", frameon=False, fontsize=9)

## Limits sized to the actual data (not the background's own long tail --
## the background runs out to ~23-29 kcal/mol, far beyond anything any of
## these 30 mutations reach, and using its 99.5th percentile here would
## waste most of the plot on empty space).
ax_scatter.set_xlim(min(-1, x_all.min()) - 0.5, x_all.max() + 0.8)
ax_scatter.set_ylim(y_all.min() - 0.8, y_all.max() + 0.8)

## Marginal histograms: the basal (background) distribution for each axis
## (same as 60_binding_vs_folding_quadrant.py).
sns.histplot(x=binding_bg, bins=80, kde=True, kde_kws=dict(bw_adjust=0.3),
             color=BG_COLOR, alpha=0.6, ax=ax_histx, line_kws=dict(lw=2))
ax_histx.tick_params(axis="x", labelbottom=False)
ax_histx.set_yticks([])
ax_histx.set_xlabel(None)
ax_histx.set_ylabel(None)
for spine in ("top", "right", "left"):
    ax_histx.spines[spine].set_visible(False)

sns.histplot(y=folding_bg, bins=80, kde=True, kde_kws=dict(bw_adjust=0.3),
             color=BG_COLOR, alpha=0.6, ax=ax_histy, line_kws=dict(lw=2))
ax_histy.tick_params(axis="y", labelleft=False)
ax_histy.set_xticks([])
ax_histy.set_xlabel(None)
ax_histy.set_ylabel(None)
for spine in ("top", "right", "bottom"):
    ax_histy.spines[spine].set_visible(False)

cbar = fig.colorbar(sc, cax=cax)
cbar.set_label("Experimental ΔV½ (mV), KCNA2 alanine scan", fontsize=11)

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/experimental_binding_vs_folding_marginals.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/experimental_binding_vs_folding_marginals.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
