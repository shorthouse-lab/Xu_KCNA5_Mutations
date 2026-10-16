#!/usr/bin/env python3
"""
Full FoldX energy decomposition for the 13 observed COSMIC pore/S5/S6
mutations, matching the column set now captured for the background
(30_pore_background_foldx.py). No new FoldX run needed -- BuildModel already
computed every term for these 13 back when 22_foldx_ddg_pore_mutations.py
was written; that script just discarded everything except the total. Reads
directly from the still-on-disk Dif_*.fxout files (same source paths as
22_foldx_ddg_pore_mutations.py, verified against that script).
"""
import csv

RESULTS_DIR = "results"
ANALYSIS = "structural"

FOLDX_COLUMNS = [
    "ddG_kcal_mol", "backbone_hbond", "sidechain_hbond", "vdw", "electro",
    "solvation_polar", "solvation_hydrophobic", "vdw_clash", "entropy_sidechain",
    "entropy_mainchain", "sloop_entropy", "mloop_entropy", "cis_bond",
    "torsional_clash", "backbone_clash", "helix_dipole", "water_bridge",
    "disulfide", "electrostatic_kon", "partial_covalent", "energy_ionisation",
    "entropy_complex",
]

## The original OneDrive-synced Dif_ files intermittently time out on read
## (cloud placeholder not materialising, a local sync issue, not a data
## problem -- different files failed on different attempts). Regenerated all
## 13 fresh into /tmp/redo_<label> instead, sidestepping the sync issue
## entirely. Verified identical total ddG to the original run for every one
## before using these (deterministic FoldX BuildModel, same base PDB +
## substitution).
MUTATIONS = [
    ("S449P", "/tmp/redo_S449P/Dif_KCNA5_protein_Repair.fxout"),
    ("A455T", "/tmp/redo_A455T/Dif_KCNA5_protein_Repair.fxout"),
    ("F464C", "/tmp/redo_F464C/Dif_KCNA5_protein_Repair.fxout"),
    ("A470T", "/tmp/redo_A470T/Dif_KCNA5_protein_Repair.fxout"),
    ("A474E", "/tmp/redo_A474E/Dif_KCNA5_protein_Repair.fxout"),
    ("V496M", "/tmp/redo_V496M/Dif_KCNA5_protein_Repair.fxout"),
    ("G497D", "/tmp/redo_G497D/Dif_KCNA5_protein_Repair.fxout"),
    ("S498L", "/tmp/redo_S498L/Dif_KCNA5_protein_Repair.fxout"),
    ("G504R", "/tmp/redo_G504R/Dif_KCNA5_protein_Repair.fxout"),
    ("V514I", "/tmp/redo_V514I/Dif_KCNA5_protein_Repair.fxout"),
    ("V514L", "/tmp/redo_V514L/Dif_KCNA5_protein_Repair.fxout"),
    ("V516I", "/tmp/redo_V516I/Dif_KCNA5_protein_Repair.fxout"),
    ("R525W", "/tmp/redo_R525W/Dif_KCNA5_protein_Repair.fxout"),
]


def read_values(path):
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    last = lines[-1].split("\t")
    return [float(x) for x in last[1:1 + len(FOLDX_COLUMNS)]]


out_path = f"{RESULTS_DIR}/tables/{ANALYSIS}/observed_foldx_breakdown.csv"
with open(out_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["label"] + FOLDX_COLUMNS)
    for label, path in MUTATIONS:
        writer.writerow([label] + read_values(path))

print(f"Saved {out_path}")
