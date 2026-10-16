#!/usr/bin/env python3
"""
Do the 19 observed T1-interface-region COSMIC mutations separate into
binding-destabilising vs. folding-destabilising effects, or is it more
continuous? Rather than force an arbitrary p90-threshold category label
onto each point (the previous version of this figure), the background
("basal") distribution for each axis is shown directly as a marginal
histogram, so the reader can judge each point's position against the full
shape of what's biophysically possible at these positions, not a single
cutoff.

x-axis: binding ddG -- the larger-magnitude of the two interface
orientations (55_fill_missing_interface_orientations.py filled in whichever
orientation the original Xu pipeline hadn't tested, so this is complete for
all 19 mutations).
y-axis: folding ddG -- FoldX BuildModel on the monomer/tetramer fold itself
(58/59: saturation screen across the same 120-216 interface region).
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from adjustText import adjust_text

RESULTS_DIR = "results"
ANALYSIS = "structural"

df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_binding_vs_folding_comparison.csv")

folding_bg = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_domain_folding_saturation_ddg_progress.csv")["ddG_kcal_mol"]
binding_bg_df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/binding_ddg_expected_vs_observed.csv")
binding_bg = binding_bg_df.loc[binding_bg_df["Group"] == "Expected", "DDG"]

x = df["max_abs_binding_ddg"]
y = df["folding_ddG"]

POINT_COLOR = "#eb6834"
BG_COLOR = "#2a78d6"
## R167 (R167C/R167H) is the most recurrently mutated T1 residue in the
## cohort (3 patients) -- it sits in the densest part of the plot (10 points
## crammed near the origin), so it was present but easy to miss/misread as
## labelling a neighbour. Given a distinct colour + larger marker so it
## stands out regardless of label crowding, not just relying on adjust_text
## to place its text legibly.
HIGHLIGHT_COLOR = "#c0392b"
is_r167 = df["Mutation"].str.startswith("R167")

fig = plt.figure(figsize=(9, 9))
gs = GridSpec(2, 2, width_ratios=(5, 1.2), height_ratios=(1.2, 5),
              wspace=0.05, hspace=0.05, figure=fig)

ax_scatter = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatter)

ax_scatter.scatter(x[~is_r167], y[~is_r167], s=90, c=POINT_COLOR, edgecolors="k", linewidths=0.8, zorder=3)
ax_scatter.scatter(x[is_r167], y[is_r167], s=130, c=HIGHLIGHT_COLOR, edgecolors="k", linewidths=1.2,
                    marker="D", zorder=4, label="R167 (most recurrent T1 mutation, n=3 patients)")
ax_scatter.legend(loc="upper left", fontsize=9, frameon=True)

texts = [ax_scatter.text(xi, yi, lab, fontsize=11 if is_r167_i else 9,
                          fontweight="bold", color=HIGHLIGHT_COLOR if is_r167_i else "black")
         for xi, yi, lab, is_r167_i in zip(x, y, df["Mutation"], is_r167)]
adjust_text(texts, ax=ax_scatter, arrowprops=dict(arrowstyle="-", color="grey", lw=0.6),
            expand=(1.5, 1.9), force_text=(0.7, 1.0))

ax_scatter.axvline(0, color="#999999", lw=1, zorder=1)
ax_scatter.axhline(0, color="#999999", lw=1, zorder=1)
ax_scatter.set_xlabel("Binding ΔΔG (kcal/mol)", fontsize=12)
ax_scatter.set_ylabel("Folding ΔΔG (kcal/mol)", fontsize=12)
ax_scatter.tick_params(axis="both", labelsize=13)

ax_scatter.set_xlim(min(-1, binding_bg.quantile(0.01)), max(x.max(), binding_bg.quantile(0.995)) + 0.5)
ax_scatter.set_ylim(min(y.min(), folding_bg.quantile(0.01)) - 0.5, max(y.max(), folding_bg.quantile(0.995)) + 0.5)

## Marginal histograms: the basal (background) distribution for each axis,
## clipped to the scatter's own view limits so the shared axes line up.
## sns.histplot's own kde=True (rather than a hand-rolled gaussian_kde
## overlay) fits and rescales the KDE per-dataset automatically -- a single
## shared bandwidth didn't track both distributions' very different
## sharpness/skew well; bw_adjust narrows Scott's default further, since
## both distributions have a sharp near-zero spike a plain default still
## oversmooths.
sns.histplot(x=binding_bg, bins=80, kde=True, kde_kws=dict(bw_adjust=0.3),
             color=BG_COLOR, alpha=0.6, ax=ax_histx, line_kws=dict(lw=2))
ax_histx.tick_params(axis="x", labelbottom=False)
ax_histx.set_yticks([])
ax_histx.set_xlabel(None)
for spine in ("top", "right", "left"):
    ax_histx.spines[spine].set_visible(False)
ax_histx.set_ylabel(None)

sns.histplot(y=folding_bg, bins=80, kde=True, kde_kws=dict(bw_adjust=0.3),
             color=BG_COLOR, alpha=0.6, ax=ax_histy, line_kws=dict(lw=2))
ax_histy.tick_params(axis="y", labelleft=False)
ax_histy.set_xticks([])
ax_histy.set_ylabel(None)
for spine in ("top", "right", "bottom"):
    ax_histy.spines[spine].set_visible(False)
ax_histy.set_xlabel(None)

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/t1_binding_vs_folding_quadrant.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/t1_binding_vs_folding_quadrant.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
