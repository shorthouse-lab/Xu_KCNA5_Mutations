#!/usr/bin/env python3
"""
FoldX BuildModel ddG for the FULL saturation-mutagenesis set (1805 = 19 non-WT
amino acids x 95 positions, 433-527), enumerated in
37_pore_full_saturation_enumeration.R. Only the 1236 rows NOT already covered
by the codon-reachable background/observed/literature runs (already_computed
column in the enumeration CSV) are actually run here -- everything else is
read from the existing progress files at merge time, not recomputed.

Same FoldX conventions as 30_pore_background_foldx.py (see that script for
the full rationale): heteromeric BuildModel on the shared repaired-WT base,
same individual_list.txt format, same parallel-worker/resumable-CSV pattern.

Unlike 30, the mutant tetramer PDB is copied out to a keep directory before
the work_dir is deleted (KCNA5_protein_Repair_1.pdb -> kept/<label>.pdb) --
this run is specifically feeding the planned D485-W472 / Y483-T477 H-bond
distance measurements next, which need the actual mutant structure, not just
its ddG. ~1236 x ~868KB is a real amount of disk (~1GB) but far less than
keeping every HOLE/CHAP intermediate would be, and it avoids re-running
FoldX a third time.
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
WORK_ROOT = "data_processed/foldx_mutants/pore_full_saturation"
KEEP_PDB_DIR = "data_processed/foldx_mutants/pore_full_saturation_pdbs"
ENUM_CSV = os.environ.get("ENUM_CSV", "results/tables/structural/pore_full_saturation_substitutions.csv")
PROGRESS_CSV = os.environ.get("PROGRESS_CSV", "results/tables/structural/pore_full_saturation_ddg_progress.csv")

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
    ## R's write.csv() writes logical columns as the strings "TRUE"/"FALSE".
    with open(ENUM_CSV, newline="") as f:
        all_rows = list(csv.DictReader(f))
    return [r for r in all_rows if r["already_computed"] == "FALSE"]


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

    try:
        result = subprocess.run(
            [FOLDX_BIN, "--command=BuildModel", "--pdb=KCNA5_protein_Repair.pdb",
             "--pdb-dir=.", "--output-dir=.", "--mutant-file=individual_list.txt"],
            cwd=work_dir, capture_output=True, text=True, timeout=300,
        )
    except subprocess.TimeoutExpired:
        ## Under heavy contention (several parallel jobs sharing 12 cores) a
        ## slow-but-not-actually-stuck FoldX run can occasionally exceed the
        ## timeout -- this used to be an unhandled exception that killed the
        ## entire multiprocessing.Pool (and every in-flight sibling worker's
        ## result with it), rather than just failing this one row. Caught
        ## here so one slow mutation costs one retry-on-resume, not the rest
        ## of the run.
        sys.stderr.write(f"[TIMEOUT] {label}: FoldX exceeded 300s, skipping (will retry on next resume)\n")
        shutil.rmtree(work_dir, ignore_errors=True)
        return row, None
    except Exception as e:
        sys.stderr.write(f"[ERROR] {label}: {e}\n")
        shutil.rmtree(work_dir, ignore_errors=True)
        return row, None

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

    mutant_pdb = os.path.join(work_dir, "KCNA5_protein_Repair_1.pdb")
    if values is not None and os.path.exists(mutant_pdb):
        os.makedirs(KEEP_PDB_DIR, exist_ok=True)
        shutil.copy(mutant_pdb, os.path.join(KEEP_PDB_DIR, f"{label}.pdb"))

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
        n_failed = 0
        with multiprocessing.Pool(N_WORKERS) as pool:
            for i, (row, values) in enumerate(pool.imap_unordered(run_one, todo), 1):
                if values is None:
                    ## Deliberately NOT written to PROGRESS_CSV -- a row only
                    ## counts as "done" (skippable on resume) once it has a
                    ## real result. Writing a None row here would make
                    ## already_done() treat a timeout/failure as complete
                    ## and it would never be retried.
                    n_failed += 1
                    continue
                writer.writerow([row["label"], row["position"], row["wt_aa"], row["mut_aa"]] + values)
                f.flush()
                if i % 10 == 0 or i == len(todo):
                    elapsed = time.time() - t_start
                    rate = elapsed / i  ## wall-clock per mutation, already reflects N_WORKERS parallelism
                    remaining = (len(todo) - i) * rate
                    ddg = values[0] if values else None
                    print(f"[{i}/{len(todo)}] {row['label']} ddG={ddg} "
                          f"({rate:.2f}s/mutation wall-clock, ~{remaining/60:.0f} min remaining)", flush=True)

    print(f"{n_failed} failed/timed-out this run (not recorded -- re-run this script to retry them).", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
