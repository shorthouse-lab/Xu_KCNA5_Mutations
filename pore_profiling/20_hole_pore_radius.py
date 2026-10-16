#!/usr/bin/env python3
"""
Per-residue pore-radius table, sourced from HOLE2 (Smart et al. 1993),
replacing the CHAP-based methods in 19_direct_pore_radius.py.

Why: CHAP's own radiusMin(s) profile disagreed sharply with its own OBJ
mesh at several axial heights (clearest at R525), traced to CHAP's
in-plane stochastic optimiser failing to converge there. Reimplementing
the same class of algorithm from scratch (straight-axis / centroid /
mesh-lookup, see 19_direct_pore_radius.py) never fully resolved the
disagreement -- confidence in all CHAP-derived numbers was undermined.
Switched to HOLE2 (github.com/osmart/hole2, conda-forge), a separate,
long-established implementation, run via Docker (linux/amd64 -- no
conda-forge build for ARM64 macOS).

HOLE's own lateral MC search has the same theoretical vulnerability CHAP
has (an unconstrained search can wander into open space with nothing to
stop it once it exits the modelled region), and running it directly on
the full T1+VSD+pore homology model (with no membrane/solvent to bound
it) reproduced exactly that failure: the search escaped up to 175 Angs
from the nearest atom, presumably by continuing straight past the S6
gate into the open cytoplasmic-side region past T1. Fix: run HOLE on
just the pore-forming module (chain A-D, residues 425-530, S5 through
the modelled end of S6) instead of the full multi-domain structure, so
there is no open route beyond the region of interest for the search to
escape into. This produced a stable, connected, non-escaping profile.

Validated at two landmarks: the selectivity filter minimum (0.754 Angs
= 0.075nm) matches CHAP's own original JSON value there almost exactly;
R525 (~6.1 Angs = 0.61nm) is wide, consistent with the VMD render
showing it near the vestibule's widest point (this was the original,
motivating discrepancy -- CHAP's profile lookup put R525 at a spurious
narrow constriction not present in the render).

Coverage note: the resolved HOLE path only spans chain-A residues
~459-527 -- it does not reach the S5 helix (433-458) at all. This is a
real structural fact, not a gap in the calculation: S5 packs laterally
against S6 rather than lining the central conduction pathway, so it was
never actually "on the pore" in the first place.
"""
import numpy as np
import pandas as pd

PDB = "data_processed/hole/KCNA5_pore_module.pdb"
SPH = "data_processed/hole/kcna5_hole.sph"   # HOLE's own 3D centerline + local radius at each sampled point

# A residue is only assigned a HOLE radius if its CA sits within this many
# Angstrom of (nearest centerline point's radius + distance to that point) --
# i.e. plausibly reachable from the CA by a real side chain. Calibrated
# against known landmarks: filter residues 479-484 (dist-radius = 2.3-5.8A,
# on-axis, should pass), R525 (2.5A, on-axis, should pass), vs. S5 helix
# residues (11-20A, off-axis/peripheral, should fail). 9.0A cleanly
# separates these groups.
CA_TO_SURFACE_MAX_A = 9.0

# Same validated pore axis used throughout this pipeline (13_chap_pore_analysis.R) --
# only used here to populate s_axial for 17_pore_s6_lollipop.py's separate,
# CHAP-JSON-based hydrophobicity lookup; not used for radius at all.
INIT_POS_NM = np.array([13.2852, -13.28204, 13.31009])
DIR_VEC = np.array([0.01344207, 0.9998711, 0.008781797])
DIR_VEC = DIR_VEC / np.linalg.norm(DIR_VEC)

# =============================================================================
# A. HOLE's sphere-centre (.sph) output: full 3D centerline points + local radius
# =============================================================================
sph = []
with open(SPH) as f:
    for line in f:
        if not line.startswith("ATOM"):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        radius_A = float(line[54:60])
        sph.append((xyz, radius_A))
sph_xyz = np.array([p[0] for p in sph])
sph_radius = np.array([p[1] for p in sph])
print(f"HOLE centerline: {len(sph)} points")

# =============================================================================
# B. Chain-A residue CA positions
# =============================================================================
cas = []
with open(PDB) as f:
    for line in f:
        if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "CA":
            resno = int(line[22:26])
            resid = line[17:20].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            cas.append((resno, resid, xyz))
print(f"Chain A residues: {len(cas)}")

# =============================================================================
# C. For each residue, nearest HOLE centerline point in true 3D space -- gated
#    by CA-to-surface plausibility (see CA_TO_SURFACE_MAX_A above) so residues
#    that are merely axially level with the pore, but laterally far from it
#    (e.g. the peripheral S5 helix), are correctly excluded rather than
#    spuriously assigned a nearby-in-one-dimension radius.
# =============================================================================
results = []
for resno, resid, xyz_A in cas:
    xyz_nm = xyz_A / 10.0
    s_axial = (xyz_nm - INIT_POS_NM) @ DIR_VEC
    d = np.linalg.norm(sph_xyz - xyz_A, axis=1)
    i = d.argmin()
    dist_to_centerline, local_radius = d[i], sph_radius[i]
    if dist_to_centerline - local_radius <= CA_TO_SURFACE_MAX_A:
        radius_nm = local_radius / 10.0
    else:
        radius_nm = np.nan
    results.append((resno, resid, radius_nm, s_axial))

res_map = pd.DataFrame(results, columns=["resno", "resid", "pore_radius_nm", "s_axial"])

# =============================================================================
# D. Domain annotation (same boundaries used throughout this pipeline)
# =============================================================================
domains = [
    ("Disordered (N-term)", 1, 120), ("Linker", 121, 127), ("T1", 128, 216),
    ("Linker", 217, 230), ("S0", 231, 240), ("Linker", 241, 246),
    ("S1", 247, 270), ("Disordered", 271, 321), ("S2", 322, 344),
    ("S3", 345, 376), ("Linker", 377, 391), ("S4", 392, 416),
    ("Linker", 417, 432), ("S5", 433, 456), ("Linker", 457, 467),
    ("Pore/selectivity filter", 468, 510), ("S6/C-term", 511, 613),
]
def classify(resno):
    for name, start, end in domains:
        if start <= resno <= end:
            return name
    return None
res_map["domain"] = res_map["resno"].apply(classify)

out_path = "results/tables/structural/chap_residue_pore_radius_map.csv"
res_map[["resno", "resid", "pore_radius_nm", "domain", "s_axial"]].to_csv(out_path, index=False)
print(f"\nSaved {out_path} ({len(res_map)} residues, {res_map.pore_radius_nm.notna().sum()} with a HOLE radius)")

print("\n--- Calibration check: selectivity filter (479-484), should read narrow ---")
print(res_map[(res_map.resno >= 479) & (res_map.resno <= 484)].to_string(index=False))
print("\n--- Calibration check: R525, should read wide (render shows widest point of vestibule) ---")
print(res_map[res_map.resno == 525].to_string(index=False))
print("\n--- S5 helix (433-456): expect NaN -- not on the resolved conduction pathway ---")
print(res_map[(res_map.resno >= 433) & (res_map.resno <= 456)].to_string(index=False))

# =============================================================================
# E. Regenerate the downstream mutated-residue ranking
# =============================================================================
cosmic = pd.read_csv("results/tables/structural/cosmic_kcna5_whole_protein_mutations.csv")
mut_positions = set(cosmic.position.unique())
res_map["is_mutated"] = res_map["resno"].isin(mut_positions)
res_map["narrowness_pctile"] = 100 * res_map["pore_radius_nm"].rank() / res_map["pore_radius_nm"].notna().sum()

mutated_map = res_map[res_map.is_mutated].sort_values("narrowness_pctile").reset_index(drop=True)
mutated_out = "results/tables/structural/chap_mutated_residues_pore_radius.csv"
mutated_map[["resno", "resid", "domain", "pore_radius_nm", "narrowness_pctile"]].to_csv(mutated_out, index=False)
print(f"\nSaved {mutated_out} ({len(mutated_map)} mutated residues; "
      f"{mutated_map.pore_radius_nm.notna().sum()} inside the resolved pathway)")
print(mutated_map[["resno", "resid", "domain", "pore_radius_nm", "narrowness_pctile"]].to_string(index=False))

# =============================================================================
# F. Narrowness-rank plot (mutated residues inside the resolved pathway only)
# =============================================================================
import matplotlib.pyplot as plt
import matplotlib.cm as cm

ranked = mutated_map[mutated_map.pore_radius_nm.notna()].copy()
ranked["resid_label"] = ranked["resid"] + ranked["resno"].astype(str)
domains_present = sorted(ranked["domain"].dropna().unique())
palette = {d: cm.tab20(i / max(len(domains_present) - 1, 1)) for i, d in enumerate(domains_present)}

fig, ax = plt.subplots(figsize=(8, 10))
order = ranked.sort_values("narrowness_pctile")
ax.barh(order["resid_label"], order["narrowness_pctile"],
        color=[palette[d] for d in order["domain"]], edgecolor="black", linewidth=0.5)
ax.set_xlabel("Narrowness percentile (0 = at the tightest constriction, 100 = widest)")
ax.set_title("Mutated KCNA5 residues ranked by proximity to the pore's tightest point\n"
              "(HOLE2, run on the isolated pore module)")
handles = [plt.Rectangle((0, 0), 1, 1, color=palette[d]) for d in domains_present]
ax.legend(handles, domains_present, title="Domain", loc="lower right", fontsize=8)
fig.tight_layout()
fig_out = "results/figures/structural/chap_mutated_residues_narrowness_rank.png"
fig.savefig(fig_out, dpi=300)
fig.savefig(fig_out.replace(".png", ".pdf"))
print(f"\nSaved {fig_out}")
