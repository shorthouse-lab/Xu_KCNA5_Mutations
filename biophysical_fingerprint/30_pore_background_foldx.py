#!/usr/bin/env python3
"""
FoldX BuildModel ddG for the 549 codon-reachable background substitutions
enumerated in 29_pore_background_enumeration.R (positions 433-527, the
structurally modelled range -- see that script's header for why 613 isn't
usable). Same base structure, same heteromeric (single chain-A subunit
mutated in the otherwise-WT repaired tetramer) setup, same individual_list.txt
convention as the original 13 observed mutations
(data_processed/foldx_mutants/pore_mutations/<label>/), verified against
those directories before this script was written.

Deliberately writes one row to the progress CSV per completed mutation,
flushed immediately, and skips any label already present there -- FoldX
BuildModel across ~549 mutations takes on the order of an hour (real,
in-session pilot timing: ~6s/mutation on this machine), long enough that
being interruptible/resumable matters more than for the original 13.

Each mutation's working directory is deleted immediately after its ddG is
parsed -- unlike the original 13 (kept for the downstream HOLE/CHAP rerun,
21_pore_radius_hydrophobicity_mutants.py), the background set was scoped to
skip HOLE/CHAP entirely (silent-failure risk at this scale, see chat), so
there's no later step that needs these structures kept around, and 549 sets
of repaired-tetramer PDBs would be >500MB of disposable intermediate output.
"""
import csv
import multiprocessing
import os
import shutil
import subprocess
import sys
import time

## Each BuildModel run is single-threaded and fully isolated to its own
## work_dir (own copy of the base PDB, own individual_list.txt, own FoldX
## process) -- no shared state between mutations, so this is embarrassingly
## parallel. 8 of this machine's 12 cores, leaving headroom rather than
## saturating every core. All file I/O for a given mutation (copy in, run,
## parse, delete) happens inside the worker; only the already-computed
## result list crosses back to the main process, which is the sole writer
## to PROGRESS_CSV -- avoids any concurrent-write risk without needing a
## lock.
N_WORKERS = 8

FOLDX_BIN = "/Users/davidshorthouse/Documents/foldx51MacC11/foldx"
BASE_PDB = "data_processed/foldx_mutants/KCNA5_protein_Repair.pdb"
WORK_ROOT = "data_processed/foldx_mutants/pore_background"
ENUM_CSV = os.environ.get("ENUM_CSV", "results/tables/structural/pore_background_substitutions.csv")
PROGRESS_CSV = os.environ.get("PROGRESS_CSV", "results/tables/structural/pore_background_ddg_progress.csv")

## Column names, in order, matching the Dif_*.fxout header exactly (sanitised
## for CSV/R use) -- captures FoldX's full energy decomposition, not just the
## total, since it costs nothing extra (BuildModel already computes every
## term; we were just discarding them at parse time). See chat: total ddG
## alone doesn't distinguish *why* a mutation is destabilising (steric clash
## vs. burying a new charge vs. losing packing) -- these components do.
FOLDX_COLUMNS = [
    "ddG_kcal_mol", "backbone_hbond", "sidechain_hbond", "vdw", "electro",
    "solvation_polar", "solvation_hydrophobic", "vdw_clash", "entropy_sidechain",
    "entropy_mainchain", "sloop_entropy", "mloop_entropy", "cis_bond",
    "torsional_clash", "backbone_clash", "helix_dipole", "water_bridge",
    "disulfide", "electrostatic_kon", "partial_covalent", "energy_ionisation",
    "entropy_complex",
]


def load_enum_rows():
    with open(ENUM_CSV, newline="") as f:
        return list(csv.DictReader(f))


def already_done():
    if not os.path.exists(PROGRESS_CSV):
        return set()
    with open(PROGRESS_CSV, newline="") as f:
        return {row["label"] for row in csv.DictReader(f)}


def run_one(row):
    label = row["label"]
    wt_aa, position, mut_aa = row["wt_aa"], row["position"], row["mut_aa"]
    work_dir = os.path.join(WORK_ROOT, label)
    os.makedirs(work_dir, exist_ok=True)
    shutil.copy(BASE_PDB, os.path.join(work_dir, "KCNA5_protein_Repair.pdb"))
    with open(os.path.join(work_dir, "individual_list.txt"), "w") as f:
        f.write(f"{wt_aa}A{position}{mut_aa};\n")

    result = subprocess.run(
        [FOLDX_BIN, "--command=BuildModel", "--pdb=KCNA5_protein_Repair.pdb",
         "--pdb-dir=.", "--output-dir=.", "--mutant-file=individual_list.txt"],
        cwd=work_dir, capture_output=True, text=True, timeout=120,
    )

    dif_path = os.path.join(work_dir, "Dif_KCNA5_protein_Repair.fxout")
    values = None
    if os.path.exists(dif_path):
        with open(dif_path) as f:
            lines = [l for l in f if l.strip()]
        if lines:
            last = lines[-1].split("\t")
            try:
                values = [float(x) for x in last[1:1 + len(FOLDX_COLUMNS)]]
            except (IndexError, ValueError):
                values = None

    shutil.rmtree(work_dir, ignore_errors=True)

    if values is None:
        sys.stderr.write(f"[FAILED] {label}: no ddG parsed. stderr tail:\n{result.stderr[-500:]}\n")
    return row, values


def main():
    rows = load_enum_rows()
    done = already_done()
    todo = [r for r in rows if r["label"] not in done]
    print(f"{len(rows)} total background substitutions, {len(done)} already done, {len(todo)} to run "
          f"({N_WORKERS} parallel workers)", flush=True)

    write_header = not os.path.exists(PROGRESS_CSV)
    with open(PROGRESS_CSV, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["label", "position", "wt_aa", "mut_aa"] + FOLDX_COLUMNS)
            f.flush()

        t_start = time.time()
        with multiprocessing.Pool(N_WORKERS) as pool:
            for i, (row, values) in enumerate(pool.imap_unordered(run_one, todo), 1):
                writer.writerow([row["label"], row["position"], row["wt_aa"], row["mut_aa"]] +
                                 (values if values is not None else [None] * len(FOLDX_COLUMNS)))
                f.flush()
                if i % 10 == 0 or i == len(todo):
                    elapsed = time.time() - t_start
                    rate = elapsed / i  ## wall-clock per mutation, already reflects N_WORKERS parallelism
                    remaining = (len(todo) - i) * rate
                    ddg = values[0] if values else None
                    print(f"[{i}/{len(todo)}] {row['label']} ddG={ddg} "
                          f"({rate:.2f}s/mutation wall-clock, ~{remaining/60:.0f} min remaining)", flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
