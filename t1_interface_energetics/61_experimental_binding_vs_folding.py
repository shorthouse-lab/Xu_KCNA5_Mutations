#!/usr/bin/env python3
"""
Extends 60_binding_vs_folding_quadrant.py's binding-vs-folding plane to the
35-mutation KCNA2 alanine-scan reference set (28_experimental_ddg_
correlation.py) rather than just the 19 observed COSMIC mutations, coloured
by the experimental activation-voltage shift (deltaV1/2) where it exists.
Motivated by chat: does a mutation's position on the binding/folding plane
relate at all to its actual functional consequence (not just whether it's
destabilising in silico), now that both a folding axis (58/59 saturation
screen) and a complete two-orientation binding axis (55) exist for this
region?

3 mutations (Q93A*, G41A*, R82A*) are expression-null/unmeasurable in the
source data (x encoded as 0, a placeholder -- see 28's docstring) --
plotted as grey open squares, not coloured as if x=0 were a real measured
zero-effect. 3 more (P105A/P192A, I108A/I195A, E111A/E198A) have no FoldX
binding ddG match at all (28's own "unmatched" set) and are dropped here
too, same as 28 dropped them from its correlation stats.
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

RESULTS_DIR = "results"
ANALYSIS = "structural"

exp = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/experimental_ddg_correlation.csv")
folding_bg = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_domain_folding_saturation_ddg_progress.csv")
binding_bg = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/binding_ddg_expected_vs_observed.csv")
binding_bg = binding_bg[binding_bg["Group"] == "Expected"].copy()
binding_bg["position"] = binding_bg["Mutation"].str.extract(r"(\d+)").astype(int)
binding_bg["mut_aa"] = binding_bg["Mutation"].str[-1]

parsed = exp["kcna5_label"].str.extract(r"^([A-Z])(\d+)([A-Z])$")
exp["position"] = parsed[1].astype(int)
exp["mut_aa"] = parsed[2]
exp["unmeasurable"] = exp["Mutation"].str.endswith("*")

## Both binding orientations, looked up fresh from the full background
## (not just the 19-mutation table 55/60 built, which only covers the
## observed COSMIC set) -- pivot so each mutation gets an A and a D column.
binding_wide = (
    binding_bg[binding_bg["position"].isin(exp["position"]) & binding_bg["mut_aa"].isin(exp["mut_aa"])]
    .pivot_table(index=["position", "mut_aa"], columns="Chain", values="DDG", aggfunc="first")
    .rename(columns={"A": "binding_AD", "D": "binding_DA"})
    .reset_index()
)
exp = exp.merge(binding_wide, on=["position", "mut_aa"], how="left")
exp = exp.merge(folding_bg[["position", "mut_aa", "ddG_kcal_mol"]].rename(columns={"ddG_kcal_mol": "folding_ddG"}),
                 on=["position", "mut_aa"], how="left")
exp["max_abs_binding_ddg"] = exp[["binding_AD", "binding_DA"]].abs().max(axis=1)

usable = exp[~exp["unmeasurable"] & exp["max_abs_binding_ddg"].notna() & exp["folding_ddG"].notna()].copy()
unmeasurable = exp[exp["unmeasurable"] & exp["max_abs_binding_ddg"].notna() & exp["folding_ddG"].notna()].copy()
unmatched = exp[exp["max_abs_binding_ddg"].isna() | exp["folding_ddG"].isna()]
print(f"Usable (has experimental x, both FoldX axes): {len(usable)}")
print(f"Unmeasurable (expression-null placeholder): {len(unmeasurable)}")
print(f"Unmatched (missing a FoldX axis): {len(unmatched)} -- {unmatched['kcna5_label'].tolist()}")

out_csv = f"{RESULTS_DIR}/tables/{ANALYSIS}/experimental_binding_vs_folding.csv"
exp.to_csv(out_csv, index=False)
print(f"Saved {out_csv}")

fig, ax = plt.subplots(figsize=(9, 7.5))

sc = ax.scatter(usable["max_abs_binding_ddg"], usable["folding_ddG"], c=usable["x"],
                 cmap="coolwarm", vmin=-15, vmax=15, s=100, edgecolors="k", linewidths=0.8, zorder=3)
ax.scatter(unmeasurable["max_abs_binding_ddg"], unmeasurable["folding_ddG"], facecolors="none",
           edgecolors="#555555", linewidths=1.5, marker="s", s=100, zorder=3,
           label="Expression-null (unmeasurable)")

cbar = fig.colorbar(sc, ax=ax)
cbar.set_label("Experimental ΔV½ (mV), KCNA2 alanine scan", fontsize=11)

texts = [ax.text(row.max_abs_binding_ddg, row.folding_ddG, row.kcna5_label, fontsize=8, fontweight="bold")
         for row in pd.concat([usable, unmeasurable]).itertuples()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.6),
            expand=(1.3, 1.6), force_text=(0.5, 0.8))

ax.axvline(0, color="#999999", lw=1, zorder=1)
ax.axhline(0, color="#999999", lw=1, zorder=1)
ax.set_xlabel("Binding ΔΔG (kcal/mol)\n(larger-magnitude of the two interface orientations)", fontsize=12)
ax.set_ylabel("Folding ΔΔG (kcal/mol)", fontsize=12)
ax.tick_params(axis="both", labelsize=12)
ax.legend(loc="lower right", frameon=False, fontsize=9)

fig.tight_layout()
out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/experimental_binding_vs_folding.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/experimental_binding_vs_folding.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
