#!/usr/bin/env python3
"""
CHAP pore-facing hydrophobicity (hydro_max_delta) for the 549 background
substitutions, piloted and validated in-session before this run:
  - G504R reproduced the known hydro_max_delta almost exactly (-0.054 here
    vs. -0.056 in mutation_effect_summary.csv).
  - CHAP's own radiusMin(s) profile, by contrast, reproduced the DOCUMENTED
    FAILURE MODE on the same pilot: for G504R (a real, HOLE-verified
    -0.008nm narrowing), CHAP's radius channel said +0.0005 -- essentially
    missing the true effect. This is the same stochastic-optimiser
    non-convergence issue that motivated switching to HOLE for radius in
    the first place (13_chap_pore_analysis.R's "Correction, 2026-07").

Because that failure mode is real and CHAP's radius output cannot be
trusted at face value, but hydrophobicity was NOT simply wrong whenever
radius disagreed (G504R again: wrong radius, right hydrophobicity) -- ruling
hydrophobicity untrustworthy across the board isn't supported by the pilot
either. Resolution: compute a per-mutation QC signal (CHAP's own radius
delta at the mutated position vs. the already-validated HOLE radius delta
for the same mutation, 33_pore_background_hole.py) and report it alongside
hydro_max_delta rather than silently trusting OR discarding CHAP's
hydrophobicity output. Large disagreement flags a mutation whose CHAP run's
underlying pathway search likely didn't converge well at that position --
worth caution, not automatic exclusion (it would have excluded G504R, the
one mutation this pipeline has the most independent evidence for).

Unlike HOLE, CHAP does NOT need the structure trimmed to the pore module --
it was run on the full tetramer for the original 13 (verified: -f and -s
were the same self-referential mutant PDB, both covering residues 120-527,
all 4 chains -- confirmed in-session that -f=reference KCNA5.pdb with a
mismatched atom count against a mutant -s fails with a GROMACS selection
error; -f=-s=<mutant>.pdb avoids this entirely and is what the original 13
runs used, per the input PDBs still on disk in data_processed/chap/).

Pilot timing: ~14.5s/mutation serial (heavier than FoldX/HOLE -- Docker +
x86_64 emulation overhead dominates), ~3.7x speedup with 8 concurrent
containers (partial scaling under emulation, still worthwhile) -- ETA ~36
min for 549 at N_WORKERS=8.
"""
import csv
import json
import multiprocessing
import os
import shutil
import subprocess
import sys
import time

import numpy as np

N_WORKERS = 8

FOLDX_BIN = "/Users/davidshorthouse/Documents/foldx51MacC11/foldx"
BASE_PDB = "data_processed/foldx_mutants/KCNA5_protein_Repair.pdb"
WORK_ROOT = "data_processed/foldx_mutants/pore_background_chap"
ENUM_CSV = os.environ.get("ENUM_CSV", "results/tables/structural/pore_background_substitutions.csv")
PROGRESS_CSV = os.environ.get("PROGRESS_CSV", "results/tables/structural/pore_background_chap_progress.csv")
HOLE_PROGRESS_CSV = "results/tables/structural/pore_background_hole_progress.csv"
WT_CHAP_JSON = "data_processed/chap/kcna5_pore_repairedWT.json"
CHAP_IMAGE = "kcna5-chap:0.9.1"

INIT_POS_NM = np.array([13.2852, -13.28204, 13.31009])
DIR_VEC = np.array([0.01344207, 0.9998711, 0.008781797])
DIR_VEC = DIR_VEC / np.linalg.norm(DIR_VEC)
HYDRO_GRID = np.linspace(-3.8, 1.8, 300)  # same resolved axial range as 24_mutation_effect_summary.py


def load_hole_deltas():
    with open(HOLE_PROGRESS_CSV, newline="") as f:
        return {row["label"]: float(row["radius_delta_nm"]) for row in csv.DictReader(f)
                if row["radius_delta_nm"] not in ("", "None")}


def ca_position(pdb_path, resno, chain="A"):
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
                if int(line[22:26]) == resno:
                    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return None


def load_enum_rows():
    with open(ENUM_CSV, newline="") as f:
        return list(csv.DictReader(f))


def already_done():
    if not os.path.exists(PROGRESS_CSV):
        return set()
    with open(PROGRESS_CSV, newline="") as f:
        return {row["label"] for row in csv.DictReader(f)}


_WT_HYDRO = None  # loaded once per worker
_WT_RADIUS = None


def run_one(row):
    global _WT_HYDRO, _WT_RADIUS
    if _WT_HYDRO is None:
        wt = json.load(open(WT_CHAP_JSON))["pathwayProfile"]
        _WT_HYDRO = np.interp(HYDRO_GRID, wt["s"], wt["pfHydrophobicityMean"])
        _WT_RADIUS = (np.array(wt["s"]), np.array(wt["radiusMin"]))

    label = row["label"]
    wt_aa, position, mut_aa = row["wt_aa"], int(row["position"]), row["mut_aa"]
    work_dir = os.path.join(WORK_ROOT, label)
    os.makedirs(work_dir, exist_ok=True)

    hydro_max_delta, chap_radius_delta = None, None
    try:
        shutil.copy(BASE_PDB, os.path.join(work_dir, "KCNA5_protein_Repair.pdb"))
        with open(os.path.join(work_dir, "individual_list.txt"), "w") as f:
            f.write(f"{wt_aa}A{position}{mut_aa};\n")

        foldx_result = subprocess.run(
            [FOLDX_BIN, "--command=BuildModel", "--pdb=KCNA5_protein_Repair.pdb",
             "--pdb-dir=.", "--output-dir=.", "--mutant-file=individual_list.txt"],
            cwd=work_dir, capture_output=True, text=True, timeout=120,
        )
        mutant_pdb = os.path.join(work_dir, "KCNA5_protein_Repair_1.pdb")
        if not os.path.exists(mutant_pdb):
            sys.stderr.write(f"[FAILED FoldX] {label}: {foldx_result.stderr[-300:]}\n")
            return row, None, None

        abs_work = os.path.abspath(work_dir)
        chap_result = subprocess.run(
            ["docker", "run", "--rm", "--platform", "linux/amd64",
             "-v", f"{abs_work}:/out", "-w", "/out", CHAP_IMAGE,
             "-f", "/out/KCNA5_protein_Repair_1.pdb", "-s", "/out/KCNA5_protein_Repair_1.pdb",
             "-sel-pathway", "Protein",
             "-pf-chan-dir-vec", *[str(x) for x in DIR_VEC],
             "-pf-init-probe-pos", *[str(x) for x in INIT_POS_NM],
             "-pf-max-free-dist", "5", "-out-filename", "chap_out"],
            capture_output=True, text=True, timeout=120,
        )
        json_path = os.path.join(work_dir, "chap_out.json")
        if not os.path.exists(json_path):
            sys.stderr.write(f"[FAILED CHAP] {label}: {chap_result.stdout[-300:]} {chap_result.stderr[-300:]}\n")
            return row, None, None

        prof = json.load(open(json_path))["pathwayProfile"]
        h = np.interp(HYDRO_GRID, prof["s"], prof["pfHydrophobicityMean"])
        diff = h - _WT_HYDRO
        hydro_max_delta = float(diff[np.abs(diff).argmax()])

        ca = ca_position(mutant_pdb, position)
        s_axial = (ca / 10.0 - INIT_POS_NM) @ DIR_VEC
        m_s, m_r = np.array(prof["s"]), np.array(prof["radiusMin"])
        chap_mut_r = np.interp(s_axial, m_s, m_r)
        chap_wt_r = np.interp(s_axial, _WT_RADIUS[0], _WT_RADIUS[1])
        chap_radius_delta = float(chap_mut_r - chap_wt_r)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return row, hydro_max_delta, chap_radius_delta


def main():
    rows = load_enum_rows()
    hole_deltas = load_hole_deltas()
    done = already_done()
    todo = [r for r in rows if r["label"] not in done]
    print(f"{len(rows)} total, {len(done)} already done, {len(todo)} to run ({N_WORKERS} workers)", flush=True)

    write_header = not os.path.exists(PROGRESS_CSV)
    with open(PROGRESS_CSV, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["label", "position", "wt_aa", "mut_aa",
                              "hydro_max_delta", "chap_radius_delta", "hole_radius_delta", "qc_disagreement"])
            f.flush()

        t_start = time.time()
        with multiprocessing.Pool(N_WORKERS) as pool:
            for i, (row, hydro, chap_r) in enumerate(pool.imap_unordered(run_one, todo), 1):
                hole_r = hole_deltas.get(row["label"])
                qc = abs(chap_r - hole_r) if (chap_r is not None and hole_r is not None) else None
                writer.writerow([row["label"], row["position"], row["wt_aa"], row["mut_aa"],
                                  hydro, chap_r, hole_r, qc])
                f.flush()
                if i % 10 == 0 or i == len(todo):
                    elapsed = time.time() - t_start
                    rate = elapsed / i
                    remaining = (len(todo) - i) * rate
                    print(f"[{i}/{len(todo)}] {row['label']} hydro_max_delta={hydro} qc_disagreement={qc} "
                          f"({rate:.2f}s/mutation wall-clock, ~{remaining/60:.0f} min remaining)", flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
