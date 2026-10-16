#!/usr/bin/env python3
"""
Standalone WT HOLE pore-radius profile, vertical orientation, meant to be
placed alongside a VMD structure render of the same channel -- extracellular
/ filter at the top, intracellular / S6 gate at the bottom, matching how a
K+ channel is conventionally rendered (membrane normal vertical, filter up).

Same axis convention as pore_radius_hydrophobicity_mutants.png (21_pore_
radius_hydrophobicity_mutants.py): s_axial = position along the validated
pore axis (13_chap_pore_analysis.R), in nm, extracellular/filter negative,
intracellular/S6 gate positive. A few residue landmarks are labelled
directly on the axis so the plot can be visually lined up against a render
even without pixel-exact cropping.
"""
import numpy as np
import matplotlib.pyplot as plt

RESULTS_DIR = "results"
ANALYSIS = "structural"
K_IONIC_RADIUS_NM = 0.15
AXIAL_YLIM = (-3.8, 1.8)  # same structurally-resolved range as 21_...py

INIT_POS_NM = np.array([13.2852, -13.28204, 13.31009])
DIR_VEC = np.array([0.01344207, 0.9998711, 0.008781797])
DIR_VEC = DIR_VEC / np.linalg.norm(DIR_VEC)


def load_sph(path):
    xyz, radius = [], []
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            radius.append(float(line[54:60]))
    return np.array(xyz), np.array(radius)


xyz, radius = load_sph("data_processed/hole/mutants/wt_hole.sph")
s_axial = (xyz / 10.0 - INIT_POS_NM) @ DIR_VEC
order = np.argsort(s_axial)
s_axial, radius_nm = s_axial[order], radius[order] / 10.0

fig, ax = plt.subplots(figsize=(3.2, 9))
ax.plot(radius_nm, s_axial, color="black", linewidth=2.6)

ax.axvline(K_IONIC_RADIUS_NM, color="firebrick", linestyle=":", linewidth=3, zorder=1)
ax.text(K_IONIC_RADIUS_NM + 0.05, AXIAL_YLIM[1] - 0.2, "K$^+$ ionic radius (0.15 nm)",
        fontsize=13, fontweight="bold", color="firebrick", va="bottom", ha="left", rotation=90,
        zorder=2, bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

ax.set_xlabel("WT HOLE pore\nradius (nm)", fontsize=16)
ax.set_xlim(1.3, 0)  # flipped: 0 on the right
ax.set_ylim(AXIAL_YLIM[1], AXIAL_YLIM[0])  # inverted: extracellular/filter at top
ax.set_yticks([])
ax.set_ylabel(None)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="both", width=1.5, labelsize=10)
for spine in ("top", "right", "bottom"):
    ax.spines[spine].set_linewidth(1.5)

fig.tight_layout()
out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/wt_pore_radius_vertical.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/wt_pore_radius_vertical.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print(f"Saved {out_png} and {out_pdf}")
