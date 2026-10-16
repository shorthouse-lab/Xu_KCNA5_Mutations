#!/usr/bin/env python3
"""
The 19 observed COSMIC T1-region mutations (60_binding_vs_folding_quadrant.py's
set), with marginal histograms showing each axis's basal/background
distribution (as in 60/63), coloured by the KCNA2 alanine-scan experimental
deltaV1/2 at the SAME RESIDUE POSITION -- broader than
62_observed_binding_vs_folding_experimental_colour.py, which required an
EXACT substitution match (only D166N qualified there). Here any of the 19
observed mutations whose position was tested at all in the KCNA2 screen gets
coloured, regardless of which substitution was tested there.

Checked directly: KCNA2 alanine-scan positions do not include 167 under any
matching -- the screen tested 166 (D166N) and 168 (N168A) but not 167, so
R167C/R167H will always be grey/uncoloured here, same as in the exact-match
version. This is a genuine data gap (position never tested), not a labelling
bug -- R167C/R167H still get a text label like every other point, just no
colour.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from adjustText import adjust_text

RESULTS_DIR = "results"
ANALYSIS = "structural"
BG_COLOR = "#2a78d6"

df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_binding_vs_folding_comparison.csv")
exp = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/experimental_ddg_correlation.csv")
folding_bg = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_domain_folding_saturation_ddg_progress.csv")["ddG_kcal_mol"]
binding_bg_df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/binding_ddg_expected_vs_observed.csv")
binding_bg = binding_bg_df.loc[binding_bg_df["Group"] == "Expected", "DDG"]

exp["kcna5_position"] = exp["kcna5_label"].str.extract(r"(\d+)").astype(int)
## position has no duplicate KCNA2 entries (checked), so a simple merge on
## position is unambiguous -- one experimental value per position at most.
pos_match = exp[["kcna5_position", "x", "kcna5_label"]].rename(
    columns={"x": "exp_value_samepos", "kcna5_label": "exp_mutation_samepos"})
df = df.merge(pos_match, left_on="position", right_on="kcna5_position", how="left")

n_matched = df["exp_value_samepos"].notna().sum()
print(f"n observed mutations: {len(df)}; n with a same-position KCNA2 match: {n_matched}")
print(df.loc[df["exp_value_samepos"].notna(), ["Mutation", "exp_mutation_samepos", "exp_value_samepos"]].to_string(index=False))
print("Positions with NO KCNA2 coverage at all:",
      df.loc[df["exp_value_samepos"].isna(), "Mutation"].tolist())

x = df["max_abs_binding_ddg"]
y = df["folding_ddG"]
has_exp = df["exp_value_samepos"].notna()

fig = plt.figure(figsize=(12, 10))
gs = GridSpec(2, 3, width_ratios=(5, 1.2, 0.3), height_ratios=(1.2, 5),
              wspace=0.08, hspace=0.05, figure=fig)

ax_scatter = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatter)
cax = fig.add_subplot(gs[1, 2])

ax_scatter.scatter(x[~has_exp], y[~has_exp], c="#cccccc", s=90, edgecolors="k", linewidths=0.8, zorder=3,
                    label="Position not tested in KCNA2 scan")
sc = ax_scatter.scatter(x[has_exp], y[has_exp], c=df.loc[has_exp, "exp_value_samepos"], cmap="coolwarm",
                         vmin=-15, vmax=15, s=140, edgecolors="k", linewidths=1.2, zorder=4,
                         label="Same-position KCNA2 match")

texts = [ax_scatter.text(xi, yi, lab, fontsize=11, fontweight="bold", ha="center", va="bottom")
         for xi, yi, lab in zip(x, y, df["Mutation"])]
adjust_text(texts, ax=ax_scatter, arrowprops=dict(arrowstyle="-", color="grey", lw=0.6),
            expand=(1.3, 1.6), force_text=(0.5, 0.8))

ax_scatter.axvline(0, color="#999999", lw=1, zorder=1)
ax_scatter.axhline(0, color="#999999", lw=1, zorder=1)
ax_scatter.set_xlabel("Binding ΔΔG (kcal/mol)", fontsize=12)
ax_scatter.set_ylabel("Folding ΔΔG (kcal/mol)", fontsize=12)
ax_scatter.tick_params(axis="both", labelsize=13)
ax_scatter.legend(loc="upper left", frameon=False, fontsize=9)

ax_scatter.set_xlim(min(-1, binding_bg.quantile(0.01)), max(x.max(), binding_bg.quantile(0.995)) + 0.5)
ax_scatter.set_ylim(min(y.min(), folding_bg.quantile(0.01)) - 0.5, max(y.max(), folding_bg.quantile(0.995)) + 0.5)

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
cbar.set_label("Experimental ΔV½ (mV), KCNA2 alanine scan\n(same residue position)", fontsize=10)

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/observed_binding_vs_folding_samepos_marginals.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/observed_binding_vs_folding_samepos_marginals.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
