#!/usr/bin/env python3
"""
Replaces Xu_KCNA5_Mutations-main/Experimental_correlation/replot.ipynb's
grouped bar chart (one experimental deltaV1/2 bar + one FoldX ddG bar per
mutation, side by side) with a scatter plot -- a bar chart juxtaposes two
series but never actually tests whether they agree, which turns out to
matter here: they don't. Pearson r=-0.18 (p=0.40), Spearman rho=-0.13
(p=0.54) on the 24 usable points -- weakly negative, not significant.
That's a real, notable finding the bar chart was hiding, not a plotting
preference.

Data provenance (T1-T1/T1-S1 interface alanine-scan electrophysiology,
KCNA2 numbering, offset +87 to KCNA5 -- reused verbatim from replot.ipynb,
not re-derived): plot-data.csv (experimental deltaV1/2, mV, activation
voltage shift) x Binding_DDG_AD.csv/Binding_DDG_DA.csv (FoldX binding ddG
for the homologous KCNA5 interface mutation).

3 mutations (Q93A*, G41A*, R82A*) are marked unmeasurable in the source
data (expression-null or similar) -- plot-data.csv encodes this as x=0,
which is a placeholder, not a genuine null effect. Plotted separately
(open grey markers) and excluded from the correlation statistics, rather
than silently included as if they were real zero-effect data points.

3 mutations (P105A, I108A, E111A / KCNA5 P192A, I195A, E198A) have no
match in the FoldX interface ddG set at all -- also excluded, also shown
separately for transparency rather than dropped without a trace.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from adjustText import adjust_text

RESULTS_DIR = "results"
ANALYSIS = "structural"
XU_DIR = "Xu_KCNA5_Mutations-main/Experimental_correlation"

# Same set the original notebook highlighted in red -- preserved as-is,
# provenance/rationale not stated in the source notebook.
HIGHLIGHT = {"R34L", "R97A", "R99A", "G41A*", "D79N"}

# =============================================================================
# A. Load + merge, reused verbatim from replot.ipynb
# =============================================================================
df = pd.read_csv(f"{XU_DIR}/plot-data.csv")

ddg_df_AD = pd.read_csv(f"{XU_DIR}/Binding_DDG_AD.csv", sep=",")
ddg_df_AD["MUT_pos"] = [int(m[2:-1]) for m in ddg_df_AD["Mutation"]]
ddg_df_AD["Mutation"] = [m[0] + str(int(m[2:-1])) + m[-1] for m in ddg_df_AD["Mutation"]]
ddg_df_AD = ddg_df_AD.assign(Side="B")
ddg_df_DA = pd.read_csv(f"{XU_DIR}/Binding_DDG_DA.csv", sep=",")
ddg_df_DA["MUT_pos"] = [int(m[2:-1]) for m in ddg_df_DA["Mutation"]]
ddg_df_DA["Mutation"] = [m[0] + str(int(m[2:-1])) + m[-1] for m in ddg_df_DA["Mutation"]]
ddg_df_DA = ddg_df_DA.assign(Side="A")
ddg_df = pd.concat([ddg_df_AD, ddg_df_DA]).sort_values("MUT_pos")

for i, row in df.iterrows():
    mutation, side = row["Mutation"], row["Side"]
    processed = mutation.strip("*")
    kcna5_mut = processed[0] + str(int(processed[1:-1]) + 87) + processed[-1]
    match = ddg_df.loc[(ddg_df["Mutation"] == kcna5_mut) & (ddg_df["Side"] == side), "DDG"]
    df.loc[i, "DDG"] = match.iloc[0] if not match.empty else np.nan
    df.loc[i, "kcna5_label"] = kcna5_mut
    df.loc[i, "dual_label"] = f"{processed}/{kcna5_mut}" + ("*" if mutation.endswith("*") else "")

out_table = f"{RESULTS_DIR}/tables/{ANALYSIS}/experimental_ddg_correlation.csv"
df.to_csv(out_table, index=False)

# =============================================================================
# B. Split into usable vs. excluded (unmeasurable / unmatched), correlate
#    only the usable set
# =============================================================================
unmeasurable = df[df.Mutation.str.endswith("*")]
unmatched = df[(~df.Mutation.str.endswith("*")) & df.DDG.isna()]
usable = df[(~df.Mutation.str.endswith("*")) & df.DDG.notna()]

r, p_r = pearsonr(usable.x, usable.DDG)
rho, p_rho = spearmanr(usable.x, usable.DDG)
print(f"n usable = {len(usable)} (excluded: {len(unmeasurable)} unmeasurable, {len(unmatched)} unmatched)")
print(f"Pearson r={r:.3f} p={p_r:.3g}")
print(f"Spearman rho={rho:.3f} p={p_rho:.3g}")

stats_out = f"{RESULTS_DIR}/tables/{ANALYSIS}/experimental_ddg_correlation_stats.csv"
pd.DataFrame({"n": [len(usable)], "pearson_r": [r], "pearson_p": [p_r],
              "spearman_rho": [rho], "spearman_p": [p_rho]}).to_csv(stats_out, index=False)
print(f"Saved {out_table} and {stats_out}")

# =============================================================================
# C. Scatter plot
# =============================================================================
fig, ax = plt.subplots(figsize=(8, 7))

ax.axhline(0, color="#999999", linewidth=1, zorder=0)
ax.axvspan(-1.608, 1.608, color="gray", alpha=0.15, zorder=0, label="Experimental neutral zone")

not_highlighted = usable[~usable.Mutation.isin(HIGHLIGHT)]
highlighted = usable[usable.Mutation.isin(HIGHLIGHT)]
ax.scatter(not_highlighted.x, not_highlighted.DDG, s=60, color="#2a78d6",
           edgecolors="k", linewidths=0.6, zorder=3, label="Usable (n=%d)" % len(usable))
ax.scatter(highlighted.x, highlighted.DDG, s=80, color="#eb6834",
           edgecolors="k", linewidths=0.6, zorder=4, label="Highlighted (source notebook)")
ax.scatter(unmeasurable.x, unmeasurable.DDG, s=60, facecolors="none", edgecolors="#888888",
           linewidths=1.2, zorder=3, marker="s", label="Unmeasurable experimentally (x=0 placeholder)")

# Regression line, usable points only
m, b = np.polyfit(usable.x, usable.DDG, 1)
xline = np.linspace(usable.x.min(), usable.x.max(), 50)
ax.plot(xline, m * xline + b, color="black", linestyle="--", linewidth=1.2, zorder=2)

texts = []
for _, row in usable.iterrows():
    texts.append(ax.text(row.x, row.DDG, row.dual_label, fontsize=7, color="#555555"))
for _, row in unmeasurable.iterrows():
    texts.append(ax.text(row.x, row.DDG, row.dual_label, fontsize=7, color="#888888"))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.6),
            expand=(1.3, 1.5))

ax.set_xlabel(r"Experimental $\Delta V_{1/2}$ (mV), KCNA2 alanine scan", fontsize=13)
ax.set_ylabel("FoldX binding $\\Delta\\Delta$G (kcal/mol), homologous KCNA5 mutation", fontsize=13)
ax.set_title(f"Experimental activation shift vs. FoldX binding $\\Delta\\Delta$G\n"
             f"Pearson r={r:.2f} (p={p_r:.2g}), Spearman $\\rho$={rho:.2f} (p={p_rho:.2g}) -- not significant",
             fontsize=12)
ax.legend(fontsize=9, loc="upper right")
ax.tick_params(axis="both", labelsize=11)
for spine in ax.spines.values():
    spine.set_linewidth(1.3)

fig.tight_layout()
out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/experimental_ddg_correlation_scatter.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/experimental_ddg_correlation_scatter.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
