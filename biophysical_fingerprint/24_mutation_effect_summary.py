#!/usr/bin/env python3
"""
One summary figure tying together every biophysical readout computed for
the 13 S5/pore/S6 mutations across 6a (HOLE radius + CHAP hydrophobicity,
mutant vs WT), 6b (FoldX ddG stability) and 6e (WT-only elastic network
model: intrinsic flexibility + distance-corrected filter coupling) --
a mutation x metric heatmap, each column ranked 0-100 across the 13
mutations so radius/energy/correlation units (which aren't comparable to
each other directly) become visually comparable, with the actual value
annotated in each cell so nothing is hidden behind the ranking.

Two columns (flexibility, filter coupling) are WT-only properties of the
position, not of the specific mutation -- included because Analysis 6e's
point is exactly that these predict which "biophysically silent"
mutations (6a/6b) might still matter through a route those methods can't
see. A third WT-only column (RASA, from Analysis 8) is included for the
same reason it explains ddG, not because it's itself a mutation effect.

Charge change and side-chain volume change need no structure at all --
computed directly from WT/mutant amino-acid identity (formal charge;
Zamyatnin 1972 residue volumes) -- added so the heatmap shows *which*
physical property plausibly drives a mutation's ddG (steric clash vs.
burying a new charge vs. losing packing), rather than leaving that as
prose in the summary doc.

Radius change uses pore_radius_nm_centerline (21_pore_radius_hydrophobicity_
mutants.py), not the CA-to-surface-gated pore_radius_nm -- the gated column
is NaN for 7 of 13 mutations (by design: most residues here aren't
pore-facing, see 20_hole_pore_radius.py), which isn't a missing-data
problem to work around, it's a real "not applicable" for the gated
column's own question. But it made most of this heatmap's radius column
blank, which isn't useful here where the question is simply "did the
profile change at this axial height", a quantity that's well-defined
everywhere along the resolved pathway regardless of how gated the
residue's own side-chain proximity is.

Rows are hierarchically clustered (scipy, Ward linkage on the ranked
0-100 values, Euclidean distance) rather than just ddG-sorted, so
mutations with similar overall biophysical *profiles* -- not just
similar ddG -- end up adjacent, with a dendrogram alongside showing the
grouping directly.
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

RESULTS_DIR = "results"
ANALYSIS = "structural"

MUTATIONS = [
    ("S449P", "S449P", 449), ("A455T", "A455T", 455), ("F464C", "F464C", 464),
    ("A470T", "A470T", 470), ("A474E", "A474E", 474), ("V496M", "V496M", 496),
    ("G497D", "G497D", 497), ("S498L", "S498L", 498), ("G504R", "g504r_het", 504),
    ("V514I", "V514I", 514), ("V514L", "V514L", 514), ("V516I", "V516I", 516),
    ("R525W", "R525W", 525),
]

# --- ddG (6b) ---------------------------------------------------------------
ddg = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/foldx_ddg_pore_mutations.csv").set_index("mutation")

# --- radius delta at own residue's axial height, mutant vs WT (6a) ----------
# Deliberately "at own position", NOT "max deviation anywhere along the
# profile" (unlike hydro_max_delta below) -- tried the latter and it
# introduces real edge noise: near the plotted range's boundaries the pore
# is already wide open (~9-10A, near bulk solvent), where independent HOLE
# Monte Carlo runs pick up small stochastic differences unrelated to the
# mutation (verified: 4 unrelated mutations showed an identical "biggest
# difference" at the same boundary point -- a tell it's search noise, not
# biology). Radius effects that are real (G504R) are local to the
# mutation's own position, so "own position" both avoids that noise and
# doesn't miss anything -- unlike hydrophobicity (next block), whose real
# G504R effect is propagated downstream and would be mostly missed by an
# own-position-only reading (+0.003 at G504 itself vs -0.056 at its peak).
pore = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/hole_chap_mutant_pore_profiles.csv")
wt_radius = pore[pore.structure == "wt"].set_index("resno")["pore_radius_nm_centerline"]
radius_delta = {}
for label, key, resno in MUTATIONS:
    row = pore[(pore.structure == key) & (pore.resno == resno)]
    r = row.pore_radius_nm_centerline.values[0] if len(row) else np.nan
    radius_delta[label] = r - wt_radius.get(resno, np.nan)

# --- max |hydrophobicity delta| anywhere in the resolved axial range (6a) ---
wt_chap = json.load(open("data_processed/chap/kcna5_pore_repairedWT.json"))["pathwayProfile"]
grid = np.linspace(-3.8, 1.8, 300)
wt_h = np.interp(grid, wt_chap["s"], wt_chap["pfHydrophobicityMean"])
hydro_delta = {}
for label, key, resno in MUTATIONS:
    chap_key = "g504r" if label == "G504R" else label
    d = json.load(open(f"data_processed/chap/kcna5_pore_{chap_key}.json"))["pathwayProfile"]
    h = np.interp(grid, d["s"], d["pfHydrophobicityMean"])
    diff = h - wt_h
    hydro_delta[label] = diff[np.abs(diff).argmax()]

# --- ENM: WT flexibility + distance-corrected filter-coupling residual (6e) -
flex = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/enm_flexibility_pore_positions.csv")
coupling = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/enm_filter_coupling_pore_positions.csv")
flex_by_resno = flex.set_index("resno")["flexibility_pctile"]
coupling_by_resno = coupling.set_index("resno")["residual_pctile"]

# --- WT burial (RASA), Analysis 8 -- already computed, not yet in this table
rasa_tbl = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/s5_pore_scaffold_cluster.csv")
rasa_by_resno = rasa_tbl.set_index("position")["RASA"]

# --- charge change and side-chain volume change, both from WT/mutant AA
# identity alone -- no structure needed, but genuinely new information: they
# let the heatmap show *which physical property* is driving a mutation's
# ddG (e.g. G504R's clash-dominated volume jump vs. A474E/G497D's new buried
# charge) rather than leaving that as prose in the summary doc.
FORMAL_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}  # else 0; H excluded (~neutral at pH 7.4)
# Zamyatnin (1972) mean amino-acid residue volumes, Angstrom^3 -- the
# standard, widely-cited scale for this kind of comparison.
VOLUME_A3 = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}

def charge_and_volume_delta(label):
    wt_aa, mut_aa = label[0], label[-1]
    d_charge = FORMAL_CHARGE.get(mut_aa, 0) - FORMAL_CHARGE.get(wt_aa, 0)
    d_volume = VOLUME_A3[mut_aa] - VOLUME_A3[wt_aa]
    return d_charge, d_volume

# =============================================================================
# Assemble
# =============================================================================
rows = []
for label, key, resno in MUTATIONS:
    d_charge, d_volume = charge_and_volume_delta(label)
    rows.append(dict(
        mutation=label,
        ddG_kcal_mol=ddg.loc[label, "ddG_kcal_mol"],
        radius_delta_nm=radius_delta[label],
        hydro_max_delta=hydro_delta[label],
        wt_flexibility_pctile=flex_by_resno.get(resno, np.nan),
        wt_filter_coupling_residual_pctile=coupling_by_resno.get(resno, np.nan),
        wt_rasa=rasa_by_resno.get(resno, np.nan),
        charge_delta=d_charge,
        volume_delta_A3=d_volume,
    ))
df = pd.DataFrame(rows).set_index("mutation")
df = df.reindex([m[0] for m in MUTATIONS])

out_path = f"{RESULTS_DIR}/tables/{ANALYSIS}/mutation_effect_summary.csv"
df.to_csv(out_path)
print(f"Saved {out_path}\n")
print(df.to_string())

# =============================================================================
# Rank each column 0-100 (comparable scale across units), magnitude-ranked
# for the two signed columns (see module docstring in the original version
# of this script for why: near-zero "no effect" shouldn't rank as extreme).
# =============================================================================
METRIC_LABELS = {
    "ddG_kcal_mol": "FoldX ddG\n(kcal/mol)",
    "radius_delta_nm": "HOLE radius\nchange (nm)",
    "hydro_max_delta": "Max CHAP\nhydrophobicity\nshift",
    "wt_flexibility_pctile": "WT ANM\nflexibility\npctile",
    "wt_filter_coupling_residual_pctile": "WT filter-coupling\nresidual pctile\n(dist-corrected)",
    "wt_rasa": "WT RASA\n(burial)",
    "charge_delta": "Charge\nchange",
    "volume_delta_A3": "Side-chain volume\nchange (A^3)",
}
# Columns needing magnitude-ranking (can be +ve or -ve; a big change either
# direction should rank as "extreme", not just a big positive one) --
# same rationale as radius_delta_nm/hydro_max_delta already had.
MAGNITUDE_RANKED = ("radius_delta_nm", "hydro_max_delta", "charge_delta", "volume_delta_A3")
rank_df = df.rank(pct=True) * 100
for col in MAGNITUDE_RANKED:
    rank_df[col] = df[col].abs().rank(pct=True) * 100

# =============================================================================
# Hierarchical clustering of rows (mutations) on their ranked biophysical
# profile -- groups mutations with a similar overall *pattern* across all 5
# metrics, not just similar ddG. Ward linkage needs no missing values; none
# remain now that radius uses the centerline column.
# =============================================================================
assert not rank_df.isna().any().any(), "unexpected NaN going into clustering"
Z = linkage(rank_df.values, method="ward", metric="euclidean")
dendro = dendrogram(Z, labels=rank_df.index.tolist(), no_plot=True)
order = dendro["leaves"]
df = df.iloc[order]
rank_df = rank_df.iloc[order]

# =============================================================================
# Hierarchical clustering of columns (metrics) the same way -- groups
# metrics whose ranked pattern across the 13 mutations tends to move
# together (e.g. if ddG and volume change usually agree), independent of
# their very different native units, since it's done on the same ranked
# 0-100 scale as the row clustering above.
# =============================================================================
Zc = linkage(rank_df.T.values, method="ward", metric="euclidean")
dendro_c = dendrogram(Zc, labels=rank_df.columns.tolist(), no_plot=True)
col_order = dendro_c["leaves"]
df = df.iloc[:, col_order]
rank_df = rank_df.iloc[:, col_order]

# Colour scale should only reflect cells with a genuinely meaningful value --
# percentile ranking stretches whatever range exists across the full 0-100
# colour scale even when every value is tiny, so a column that's mostly
# noise (radius_delta_nm, hydro_max_delta -- see 20_hole_pore_radius.py /
# module docstring above) still shows "strong" colours for differences that
# aren't real. Mask cells below a genuine magnitude threshold per column,
# not just a display-rounding coincidence. Thresholds chosen from the data's
# own gap, not arbitrarily: radius_delta_nm has a ~40x jump from the next-
# largest value (F464C, 2.0e-4 nm) up to G504R (8.0e-3 nm); hydro_max_delta
# has 11 of 13 mutations tied at exactly the same floor value (0.000635 --
# itself a tell that it's a shared numerical artifact, not 11 independent
# real measurements) before jumping to V516I (0.017) and G504R (0.056).
# Only these two delta columns need this -- ddG and the two WT ANM
# percentile columns are properties that genuinely vary across all 13
# mutations, not differences clustered at a shared near-zero floor.
# Clustering above already used the real (unmasked) ranks -- this is a
# display-only mask, not a data change.
LOW_VALUE_THRESHOLD = {"radius_delta_nm": 0.001, "hydro_max_delta": 0.005}
color_df = rank_df.copy()
for col, threshold in LOW_VALUE_THRESHOLD.items():
    color_df.loc[df[col].abs() < threshold, col] = np.nan

# =============================================================================
# Figure: dendrogram + heatmap side by side, same row order. Explicit
# figure-fraction axes positions rather than subplots/gridspec+tight_layout
# -- both fought the dendrogram's leaf-label placement (its rightmost leaf
# stems sit exactly at the dendrogram/heatmap boundary, colliding with
# whichever axis's tick labels ended up nearest that seam). Explicit
# coordinates give a real, guaranteed-empty corridor for the row labels.
# =============================================================================
fig = plt.figure(figsize=(15, 8.5))
ax_dendro     = fig.add_axes([0.05, 0.10, 0.12, 0.60])
ax            = fig.add_axes([0.24, 0.10, 0.66, 0.60])
ax_dendro_col = fig.add_axes([0.24, 0.71, 0.66, 0.10])
cax           = fig.add_axes([0.93, 0.10, 0.015, 0.60])

dendrogram(Zc, labels=rank_df.columns.tolist(), orientation="top", ax=ax_dendro_col,
           color_threshold=0, above_threshold_color="#555555", no_labels=True)
ax_dendro_col.set_xticks([])
ax_dendro_col.set_yticks([])
for spine in ax_dendro_col.spines.values():
    spine.set_visible(False)

dendrogram(Z, no_labels=True, orientation="left", ax=ax_dendro,
           color_threshold=0, above_threshold_color="#555555")
ax_dendro.set_xticks([])
ax_dendro.set_yticks([])
for spine in ax_dendro.spines.values():
    spine.set_visible(False)

cmap = plt.get_cmap("RdYlBu_r").copy()
cmap.set_bad(color="#dddddd")
im = ax.imshow(np.ma.masked_invalid(color_df.values), cmap=cmap, vmin=0, vmax=100, aspect="auto")

DISPLAY_FMT = {
    "radius_delta_nm": "{:+.3f}", "hydro_max_delta": "{:+.3f}",
    "charge_delta": "{:+.0f}", "volume_delta_A3": "{:+.1f}", "wt_rasa": "{:.2f}",
}
for i in range(rank_df.shape[0]):
    for j, col in enumerate(df.columns):
        val = df.iloc[i, j]
        if pd.isna(val):
            fmt, txt_colour = "n/a", "#888888"
        else:
            fmt = DISPLAY_FMT.get(col, "{:+.1f}").format(val)
            r = color_df.iloc[i, j]
            txt_colour = "white" if (not pd.isna(r) and (r > 75 or r < 25)) else "black"
        ax.text(j, i, fmt, ha="center", va="center", fontsize=9, color=txt_colour)

ax.set_xticks(range(len(df.columns)))
ax.set_xticklabels([METRIC_LABELS[c] for c in df.columns], fontsize=9)
ax.set_yticks(range(len(df.index)))
ax.set_yticklabels(df.index, fontsize=11)
fig.suptitle("Mutational effect summary, 13 S5/pore/S6 mutations\n"
             "colour = rank within column (0-100), below-threshold cells left uncoloured;\n"
             "annotation = actual value\n"
             "rows hierarchically clustered (Ward) on their ranked profile",
             fontsize=12)
cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Rank within column (percentile)", fontsize=9)

fig_out = f"{RESULTS_DIR}/figures/{ANALYSIS}/mutation_effect_summary_heatmap.png"
fig.savefig(fig_out, dpi=300)
fig.savefig(fig_out.replace(".png", ".pdf"))
print(f"\nSaved {fig_out}")
