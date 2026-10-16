#!/usr/bin/env python3
"""
Per-residue pore-radius table, sourced directly from the CHAP OBJ mesh
(kcna5_pore_vmd.obj) rather than from CHAP's own radiusMin(s) JSON profile
or from an independent reimplementation.

History, briefly (full detail in conversation, not repeated here): the
original approach (13_chap_pore_analysis.R) projected each residue onto a
straight pore axis and looked up CHAP's radiusMin(s) profile there. That
was found to disagree sharply with CHAP's own OBJ mesh at several axial
heights (clearest at R525: profile said 0.076nm, mesh said ~0.85-0.9nm) --
a structure/file mismatch was directly ruled out, pointing to CHAP's
in-plane radius optimiser (stochastic simulated annealing, per its own
documented parameters) failing to converge at those heights. Two
from-scratch reimplementations were then tried (straight-axis nearest-atom
search; 4-fold-symmetric-centroid nearest-atom search) to sidestep CHAP's
profile entirely -- both are *independent, per-height* calculations,
whereas CHAP's true centerline is the output of a *sequential,
path-dependent* search (each axial step's lateral position is optimised
starting from the previous step's converged position, not computed fresh
in isolation). That structural difference meant neither reimplementation
matched CHAP's own mesh well, especially in bent/irregular regions (the
S5-pore linker, the post-filter loop) -- correlation ~0.6-0.8, not the
near-1:1 agreement a correct method should give against ground truth.

Resolution: stop trying to reverse-engineer CHAP's undocumented sequential
optimiser, and instead just read the radius directly off CHAP's own mesh
-- the same mesh the VMD render displays. This guarantees the plot and the
render agree by construction, since they're now built from the identical
underlying geometry, and sidesteps the profile-vs-mesh discrepancy
entirely (the profile is no longer used at all).

Method: for each chain-A residue's CA, find the nearest faces (by
face-centroid distance) in the mesh's "avg_radius" material group, and
read off their material step (0-99, a linear gradient CHAP itself built
over [0, out-cap] for that CHAP run -- kcna5_pore_vmd.obj was generated
with -pf-max-free-dist 2, so step 99 = 2.0nm). Verified this calibration
against two independent checks: the filter (479-484) reads ~0.04-0.09nm,
consistent with real experimental HOLE profiles of K+ channel filters;
R525 reads ~0.85-0.9nm, consistent with a from-scratch raw-atom
nearest-neighbour calculation done independently of the mesh entirely.
"""
import numpy as np
import pandas as pd

PDB = "data_processed/chap/kcna5_pore_vmd.pdb"     # same structure the VMD render loads
OBJ = "data_processed/chap/kcna5_pore_vmd.obj"     # same mesh the VMD render displays
MAX_RADIUS_NM = 2.0    # this run's -pf-max-free-dist cap; material step 99 = this value
N_NEAREST_FACES = 10   # average over a few nearest faces for stability

# Same validated pore axis used throughout this pipeline (13_chap_pore_analysis.R) --
# only used here to keep the s_axial column downstream consumers (17_pore_s6_lollipop.py,
# for its separate, unaffected hydrophobicity lookup) rely on; not used for radius at all.
INIT_POS_NM = np.array([13.2852, -13.28204, 13.31009])
DIR_VEC = np.array([0.01344207, 0.9998711, 0.008781797])
DIR_VEC = DIR_VEC / np.linalg.norm(DIR_VEC)

# =============================================================================
# A. Chain-A residue CA positions
# =============================================================================
cas = []
with open(PDB) as f:
    for line in f:
        if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "CA":
            resno = int(line[22:26])
            resid = line[17:20].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) / 10.0
            cas.append((resno, resid, xyz))
print(f"Chain A residues: {len(cas)}")

# =============================================================================
# B. Parse the OBJ mesh: faces + material step, "avg_radius" group only
#    (confirmed identical to the plain "radius" group in this static run)
# =============================================================================
verts = []
faces = []
current_group = None
current_mtl_step = None
with open(OBJ) as f:
    for line in f:
        if line.startswith("v "):
            p = line.split()
            verts.append((float(p[1]), float(p[2]), float(p[3])))
        elif line.startswith("g "):
            current_group = line.split()[1]
        elif line.startswith("usemtl "):
            mtl = line.split()[1]
            current_mtl_step = int(mtl.split("_")[-1]) if (current_group == "avg_radius" and mtl.startswith("avg_radius_")) else None
        elif line.startswith("f ") and current_group == "avg_radius" and current_mtl_step is not None:
            idx = [int(tok.split("/")[0]) - 1 for tok in line.split()[1:4]]
            faces.append((idx[0], idx[1], idx[2], current_mtl_step))
verts = np.array(verts) / 10.0  # Angstrom -> nm
face_centroids = np.array([(verts[a] + verts[b] + verts[c]) / 3 for a, b, c, _ in faces])
face_steps = np.array([m for _, _, _, m in faces])
print(f"avg_radius faces: {len(faces)}")

# =============================================================================
# C. For each residue, nearest mesh faces -> radius
# =============================================================================
results = []
for resno, resid, ca in cas:
    d = np.linalg.norm(face_centroids - ca, axis=1)
    nearest = np.argsort(d)[:N_NEAREST_FACES]
    mean_step = face_steps[nearest].mean()
    radius_nm = mean_step / 99 * MAX_RADIUS_NM
    s_axial = (ca - INIT_POS_NM) @ DIR_VEC
    results.append((resno, resid, radius_nm, d[nearest].mean(), s_axial))

res_map = pd.DataFrame(results, columns=["resno", "resid", "pore_radius_nm", "mean_dist_to_nearest_faces_nm", "s_axial"])

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
print(f"\nSaved {out_path} ({len(res_map)} residues)")

print("\n--- Calibration check: selectivity filter (479-484), should read narrow ---")
print(res_map[(res_map.resno >= 479) & (res_map.resno <= 484)].to_string(index=False))
print("\n--- Calibration check: R525, should read wide (~0.85-0.9nm) ---")
print(res_map[res_map.resno == 525].to_string(index=False))
print("\n--- Previously-disputed regions: S5-pore linker (452-469), post-filter loop (485-494) ---")
print(res_map[((res_map.resno >= 452) & (res_map.resno <= 469)) | ((res_map.resno >= 485) & (res_map.resno <= 494))].to_string(index=False))

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
print(f"\nSaved {mutated_out} ({len(mutated_map)} mutated residues)")
print(mutated_map[["resno", "resid", "domain", "pore_radius_nm", "narrowness_pctile"]].to_string(index=False))

# =============================================================================
# F. Narrowness-rank plot
# =============================================================================
import matplotlib.pyplot as plt
import matplotlib.cm as cm

mutated_map["resid_label"] = mutated_map["resid"] + mutated_map["resno"].astype(str)
domains_present = sorted(mutated_map["domain"].dropna().unique())
palette = {d: cm.tab20(i / max(len(domains_present) - 1, 1)) for i, d in enumerate(domains_present)}

fig, ax = plt.subplots(figsize=(8, 10))
order = mutated_map.sort_values("narrowness_pctile")
ax.barh(order["resid_label"], order["narrowness_pctile"],
        color=[palette[d] for d in order["domain"]], edgecolor="black", linewidth=0.5)
ax.set_xlabel("Narrowness percentile (0 = at the tightest constriction, 100 = widest)")
ax.set_title("Mutated KCNA5 residues ranked by proximity to the pore's tightest point\n(sourced directly from the CHAP mesh, same geometry as the VMD render)")
handles = [plt.Rectangle((0, 0), 1, 1, color=palette[d]) for d in domains_present]
ax.legend(handles, domains_present, title="Domain", loc="lower right", fontsize=8)
fig.tight_layout()
fig_out = "results/figures/structural/chap_mutated_residues_narrowness_rank.png"
fig.savefig(fig_out, dpi=300)
fig.savefig(fig_out.replace(".png", ".pdf"))
print(f"\nSaved {fig_out}")
