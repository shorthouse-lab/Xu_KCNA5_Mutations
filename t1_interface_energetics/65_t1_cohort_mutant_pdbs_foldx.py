#!/usr/bin/env python3
"""
Regenerate FoldX BuildModel mutant tetramer PDBs for the 6 cohort mutations
used in the ANM elastic-network comparison (66_t1_anm_elastic_network.py):
D166N, R167C, R167H, G128R, R121H, R121C.

59_t1_domain_folding_saturation_foldx.py already ran BuildModel for all of
these (they're part of the 1691-mutation saturation screen) but deletes each
mutant work_dir immediately after parsing the ddG value -- so the ddG
numbers survive in the progress CSV but no structures do. This script is
the same BuildModel call, just for these 6 labels only, keeping the
resulting "_1.pdb" mutant model instead of discarding it, so ANM has an
actual perturbed structure (not just a WT structure with a different
sequence label) to compute on.
"""
import os
import shutil
import subprocess

FOLDX_BIN = "/Users/davidshorthouse/Documents/foldx51MacC11/foldx"
BASE_PDB = "data_processed/foldx_mutants/KCNA5_protein_Repair.pdb"
WORK_ROOT = "data_processed/foldx_mutants/t1_cohort_anm_pdbs"
OUT_DIR = "data_processed/foldx_mutants/t1_cohort_anm_pdbs_final"

## label: (wt_aa, position, mut_aa) -- chain A only, same convention as the
## saturation screen (single subunit mutated in the otherwise-WT tetramer).
MUTATIONS = {
    "D166N": ("D", 166, "N"),
    "R167C": ("R", 167, "C"),
    "R167H": ("R", 167, "H"),
    "G128R": ("G", 128, "R"),
    "R121H": ("R", 121, "H"),
    "R121C": ("R", 121, "C"),
}

os.makedirs(OUT_DIR, exist_ok=True)

for label, (wt_aa, position, mut_aa) in MUTATIONS.items():
    work_dir = os.path.join(WORK_ROOT, label)
    os.makedirs(work_dir, exist_ok=True)
    shutil.copy(BASE_PDB, os.path.join(work_dir, "KCNA5_protein_Repair.pdb"))
    with open(os.path.join(work_dir, "individual_list.txt"), "w") as f:
        f.write(f"{wt_aa}A{position}{mut_aa};\n")

    print(f"Running FoldX BuildModel for {label}...", flush=True)
    result = subprocess.run(
        [FOLDX_BIN, "--command=BuildModel", "--pdb=KCNA5_protein_Repair.pdb",
         "--pdb-dir=.", "--output-dir=.", "--mutant-file=individual_list.txt"],
        cwd=work_dir, capture_output=True, text=True, timeout=300,
    )
    mutant_pdb = os.path.join(work_dir, "KCNA5_protein_Repair_1.pdb")
    if os.path.exists(mutant_pdb):
        shutil.copy(mutant_pdb, os.path.join(OUT_DIR, f"{label}.pdb"))
        print(f"  Saved {OUT_DIR}/{label}.pdb", flush=True)
    else:
        print(f"  [FAILED] {label}: no mutant PDB produced. stderr tail:\n{result.stderr[-500:]}", flush=True)
    shutil.rmtree(work_dir, ignore_errors=True)

## Keep one copy of the WT base structure alongside the mutants for ANM to
## load from the same directory.
shutil.copy(BASE_PDB, os.path.join(OUT_DIR, "WT.pdb"))
print("Done.", flush=True)
