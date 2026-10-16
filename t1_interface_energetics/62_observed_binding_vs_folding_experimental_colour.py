#!/usr/bin/env python3
"""
61_experimental_binding_vs_folding.py restricted to just the 19 observed
COSMIC T1-region mutations (60_binding_vs_folding_quadrant.py's set),
coloured by experimental deltaV1/2 where an EXACT substitution match exists
in the KCNA2 alanine-scan reference set (28_experimental_ddg_correlation.py)
-- not just the same residue with a different substitution. Checked: of the
19 observed mutations, 7 share a residue with something in the experimental
set (121, 128, 130, 166, 184, 186, 195), but only D166N is an exact
substitution match (KCNA2 D79N). Colouring by a same-residue-different-
substitution value would silently overstate what's actually comparable, so
every other point stays grey/uncoloured.
"""
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

RESULTS_DIR = "results"
ANALYSIS = "structural"

df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_binding_vs_folding_comparison.csv")
exp = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/experimental_ddg_correlation.csv")

df = df.merge(exp[["kcna5_label", "x"]].rename(columns={"kcna5_label": "Mutation"}), on="Mutation", how="left")
print(f"n observed mutations: {len(df)}; n with an exact experimental match: {df['x'].notna().sum()}")
print(df.loc[df["x"].notna(), ["Mutation", "x"]].to_string(index=False))

x = df["max_abs_binding_ddg"]
y = df["folding_ddG"]
has_exp = df["x"].notna()

fig, ax = plt.subplots(figsize=(9, 7.5))

ax.scatter(x[~has_exp], y[~has_exp], c="#cccccc", s=90, edgecolors="k", linewidths=0.8, zorder=3,
           label="No exact experimental match")
sc = ax.scatter(x[has_exp], y[has_exp], c=df.loc[has_exp, "x"], cmap="coolwarm", vmin=-15, vmax=15,
                 s=140, edgecolors="k", linewidths=1.2, zorder=4, label="Exact experimental match")

cbar = fig.colorbar(sc, ax=ax)
cbar.set_label("Experimental ΔV½ (mV), KCNA2 alanine scan", fontsize=11)

texts = [ax.text(xi, yi, lab, fontsize=9, fontweight="bold") for xi, yi, lab in zip(x, y, df["Mutation"])]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.6),
            expand=(1.3, 1.6), force_text=(0.5, 0.8))

ax.axvline(0, color="#999999", lw=1, zorder=1)
ax.axhline(0, color="#999999", lw=1, zorder=1)
ax.set_xlabel("Binding ΔΔG (kcal/mol)\n(larger-magnitude of the two interface orientations)", fontsize=12)
ax.set_ylabel("Folding ΔΔG (kcal/mol)", fontsize=12)
ax.tick_params(axis="both", labelsize=12)
ax.legend(loc="upper left", frameon=False, fontsize=9)

fig.tight_layout()
out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/observed_binding_vs_folding_experimental_colour.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/observed_binding_vs_folding_experimental_colour.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
