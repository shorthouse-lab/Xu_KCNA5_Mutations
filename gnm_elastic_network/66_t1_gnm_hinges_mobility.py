#!/usr/bin/env python3
"""
GNM (Gaussian Network Model) on the full WT tetramer (data_processed/
foldx_mutants/t1_cohort_anm_pdbs_final/WT.pdb, chains A-D, residues
120-527) -- NOT a WT-vs-mutant comparison (65_t1_cohort_mutant_pdbs_foldx.py
already showed FoldX BuildModel leaves the Calpha trace numerically
identical to WT for all 6 cohort mutations, so an elastic network built on
those mutant structures would be a null result by construction, same
topology in, same modes out).

Instead: ask where the cohort's recurrently-mutated T1 positions (R121,
G128, D166, R167) sit relative to the WT structure's own predicted slow
global modes -- specifically, whether any of them sit at or near a HINGE
(a near-zero-fluctuation pivot point in the dominant slow mode, these are
classically where a point mutation has outsized functional leverage despite
low local flexibility) versus in a high-mobility segment that itself
contributes to the motion. This is independent, structure-only evidence for
(or against) a dynamics-mediated mechanism, complementing the D166N MD
finding (decoupled inter-chain distance correlations) without needing a
mutant-structure ANM run that we've shown can't carry any signal here.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.linalg
import prody

prody.confProDy(verbosity="none")

## ProDy 2.4.1's eigtools still calls scipy.linalg.eigh(..., turbo=, eigvals=),
## both removed in scipy>=1.14 (eigvals -> subset_by_index, turbo dropped
## entirely). Shim rather than patching the installed package.
_orig_eigh = scipy.linalg.eigh
def _eigh_compat(a, turbo=None, eigvals=None, **kwargs):
    if eigvals is not None:
        kwargs["subset_by_index"] = eigvals
    return _orig_eigh(a, **kwargs)
scipy.linalg.eigh = _eigh_compat

RESULTS_DIR = "results"
ANALYSIS = "structural"
PDB_PATH = "data_processed/foldx_mutants/t1_cohort_anm_pdbs_final/WT.pdb"

RESIDUES_OF_INTEREST = {121: "R121", 128: "G128", 166: "D166", 167: "R167"}
## T1_START moved from 128 to 120 for this and the companion lollipop figure
## (57_t1_interface_lollipop.py) only -- see that script's header for the
## rationale (UniProt's literature T1 annotation, disordered-region boundary).
## Brings R121 inside the shaded T1 span rather than sitting just left of it.
T1_START, T1_END = 120, 216

wt = prody.parsePDB(PDB_PATH)
ca = wt.select("calpha")
print(f"Tetramer: {ca.numAtoms()} Calpha atoms across chains {sorted(set(ca.getChids()))}")

gnm = prody.GNM("KCNA5_tetramer")
gnm.buildKirchhoff(ca, cutoff=10.0)
gnm.calcModes(n_modes=20)

N_SLOW_MODES = 5  # mode 1 alone turned out to be dominated by motion far
                   # outside T1 (hinges at 225, 526/527) -- checked modes
                   # 1-5 individually and found mode 2 and mode 4 instead
                   # place hinges right next to the cohort positions
                   # (mode 2: 129/164/171; mode 4: 161-163/170), so hinge
                   # proximity is assessed across this slow-mode set, not
                   # mode 1 in isolation.
sqf_all20 = prody.calcSqFlucts(gnm[:20])

chids = ca.getChids()
resnums = ca.getResnums()
df = pd.DataFrame({"chain": chids, "resnum": resnums, "sqf_20modes": sqf_all20})
df["mobility_percentile"] = df["sqf_20modes"].rank(pct=True) * 100

hinges_per_mode = {}
for m in range(N_SLOW_MODES):
    hinge_idx = prody.calcHinges(gnm[m], atoms=ca)
    hinge_resnums_A = sorted(r for c, r in zip(hinge_idx.getChids(), hinge_idx.getResnums()) if c == "A")
    hinges_per_mode[m + 1] = hinge_resnums_A
    print(f"Mode {m + 1}: chain-A hinge sites = {hinge_resnums_A}")

chainA = df[df["chain"] == "A"].reset_index(drop=True)

print("\n=== Cohort positions (chain A) relative to slow-mode hinges & mobility ===")
rows = []
for pos, label in RESIDUES_OF_INTEREST.items():
    row = chainA[chainA["resnum"] == pos].iloc[0]
    best_mode, best_hinge, best_dist = None, None, 10**9
    for mode, hinge_list in hinges_per_mode.items():
        if not hinge_list:
            continue
        dists = [abs(pos - h) for h in hinge_list]
        i = int(np.argmin(dists))
        if dists[i] < best_dist:
            best_dist, best_hinge, best_mode = dists[i], hinge_list[i], mode
    rows.append({
        "position": pos, "label": label,
        "mobility_percentile_20modes": round(row["mobility_percentile"], 1),
        "nearest_hinge_mode": best_mode, "nearest_hinge_resnum": best_hinge,
        "seq_distance_to_nearest_hinge": best_dist,
    })
summary = pd.DataFrame(rows)
print(summary.to_string(index=False))
summary.to_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/t1_cohort_positions_vs_gnm_hinges.csv", index=False)

# =============================================================================
# Figure: chain-A mobility profile (20-mode sqf) with T1 domain box,
# mode-2 hinge sites (the slow mode whose hinges sit closest to the cohort
# positions), and cohort positions marked.
# =============================================================================
HIGHLIGHT_MODE = 2
highlight_hinges = hinges_per_mode[HIGHLIGHT_MODE]

fig, ax = plt.subplots(figsize=(11, 5))
ax.plot(chainA["resnum"], chainA["sqf_20modes"], color="#2a78d6", lw=1.8, zorder=3)
ax.fill_between(chainA["resnum"], chainA["sqf_20modes"], color="#2a78d6", alpha=0.15, zorder=2)

ax.axvspan(T1_START, T1_END, color="#FF8D00", alpha=0.15, zorder=1, label="T1 domain")

for i, h in enumerate(highlight_hinges):
    ax.axvline(h, color="#555555", lw=1, linestyle="--", zorder=2,
               label=f"Mode {HIGHLIGHT_MODE} hinge site" if i == 0 else None)

## D166/R167 sit 1 residue apart, close enough that their labels collide
## directly above the points at this width -- anchor them left/right with a
## wider offset instead of centering both (same approach as
## 57_t1_interface_lollipop.py's V195M/F196L fix).
LABEL_STYLE = {"D166": dict(xytext=(-25, 12), ha="right"), "R167": dict(xytext=(25, 12), ha="left")}
for pos, label in RESIDUES_OF_INTEREST.items():
    yi = chainA.loc[chainA["resnum"] == pos, "sqf_20modes"].values[0]
    ax.scatter([pos], [yi], s=120, c="#c0392b", edgecolors="k", linewidths=1.2, zorder=5)
    style = LABEL_STYLE.get(label, dict(xytext=(0, 12), ha="center"))
    ax.annotate(label, (pos, yi), textcoords="offset points", fontsize=10, fontweight="bold",
                color="#c0392b", **style)

ax.set_xlim(T1_START - 15, T1_END + 15)
ax.set_xlabel("KCNA5 residue (chain A)", fontsize=12)
ax.set_ylabel("Mean-square fluctuation\n(GNM, slowest 20 modes, a.u.)", fontsize=12)
ax.legend(loc="upper right", fontsize=9, frameon=True)
ax.tick_params(axis="both", labelsize=11)
sns.despine(ax=ax)
fig.tight_layout()

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/t1_gnm_mobility_hinges.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/t1_gnm_mobility_hinges.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nSaved {out_png} and {out_pdf}")
