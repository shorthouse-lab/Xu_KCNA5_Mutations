#!/usr/bin/env python3
"""
Regenerates the FoldX mutant tetramer PDB for the 569 mutations whose ddG
etc. we already have (13 observed + 549 background + 7 literature-only --
2 literature labels, H463R/V476A, are already background members) but whose
structures were deleted at the time (30/34's disposable-work-dir design,
before the H-bond distance idea existed). Doesn't re-parse or store ddG --
that's already in pore_background_ddg_progress.csv / mutation_effect_summary.csv
/ literature_ctype_mutations.csv. Only the PDB is needed here, for
40_hbond_distances.py to read.

Input list built ad hoc in chat (not a committed enumeration script, since
it's just the union of three already-existing tables) -- see
/tmp/regen_pdb_list.csv at the time this was run, or regenerate the same
union directly from mutation_effect_summary.csv + pore_background_ddg_
progress.csv + literature_ctype_mutations.csv labels if re-running later.
"""
import csv
import multiprocessing
import os
import shutil
import subprocess
import sys

N_WORKERS = 4  # reduced from 8: script 38 is still running its own 8 workers concurrently; 12 total keeps this machine's 12 cores from being badly oversubscribed
FOLDX_BIN = "/Users/davidshorthouse/Documents/foldx51MacC11/foldx"
BASE_PDB = "data_processed/foldx_mutants/KCNA5_protein_Repair.pdb"
WORK_ROOT = "data_processed/foldx_mutants/regen_for_hbonds"
KEEP_PDB_DIR = "data_processed/foldx_mutants/pore_full_saturation_pdbs"  # same dir as 38 -- one combined PDB store
INPUT_CSV = os.environ.get("INPUT_CSV", "/tmp/regen_pdb_list.csv")


def run_one(row):
    label, position, wt_aa, mut_aa = row["label"], row["position"], row["wt_aa"], row["mut_aa"]
    out_pdb = os.path.join(KEEP_PDB_DIR, f"{label}.pdb")
    if os.path.exists(out_pdb):
        return label, True

    work_dir = os.path.join(WORK_ROOT, label)
    os.makedirs(work_dir, exist_ok=True)
    shutil.copy(BASE_PDB, os.path.join(work_dir, "KCNA5_protein_Repair.pdb"))
    with open(os.path.join(work_dir, "individual_list.txt"), "w") as f:
        f.write(f"{wt_aa}A{position}{mut_aa};\n")

    try:
        result = subprocess.run(
            [FOLDX_BIN, "--command=BuildModel", "--pdb=KCNA5_protein_Repair.pdb",
             "--pdb-dir=.", "--output-dir=.", "--mutant-file=individual_list.txt"],
            cwd=work_dir, capture_output=True, text=True, timeout=300,
        )
    except subprocess.TimeoutExpired:
        ## Same fix as 38_pore_full_saturation_foldx.py: an uncaught timeout
        ## here used to kill the whole Pool. This script's resume check
        ## (out_pdb existence, above) is unaffected by a mid-run crash
        ## either way, but catching it means one slow mutation doesn't take
        ## the rest of the batch down with it.
        sys.stderr.write(f"[TIMEOUT] {label}\n")
        shutil.rmtree(work_dir, ignore_errors=True)
        return label, False

    mutant_pdb = os.path.join(work_dir, "KCNA5_protein_Repair_1.pdb")
    ok = os.path.exists(mutant_pdb)
    if ok:
        os.makedirs(KEEP_PDB_DIR, exist_ok=True)
        shutil.copy(mutant_pdb, out_pdb)
    else:
        sys.stderr.write(f"[FAILED] {label}: {result.stderr[-300:]}\n")
    shutil.rmtree(work_dir, ignore_errors=True)
    return label, ok


def main():
    with open(INPUT_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    print(f"{len(rows)} PDBs to regenerate ({N_WORKERS} workers)", flush=True)
    n_ok = 0
    with multiprocessing.Pool(N_WORKERS) as pool:
        for i, (label, ok) in enumerate(pool.imap_unordered(run_one, rows), 1):
            n_ok += ok
            if i % 50 == 0 or i == len(rows):
                print(f"[{i}/{len(rows)}] ok so far: {n_ok}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
