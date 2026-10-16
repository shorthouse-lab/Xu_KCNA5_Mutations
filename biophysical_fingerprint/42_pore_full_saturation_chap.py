#!/usr/bin/env python3
"""
CHAP pore-facing hydrophobicity for the full saturation set -- extends
34_pore_background_chap.py's validated method (self-referential -f/-s on the
full mutant tetramer, max |hydrophobicity delta| vs WT across the resolved
axial range, plus the CHAP-radius-vs-HOLE-radius QC disagreement signal)
from the 549 codon-reachable background to every possible substitution.

Like 41_pore_full_saturation_hole.py, reuses the already-saved mutant PDBs
(data_processed/foldx_mutants/pore_full_saturation_pdbs/) instead of
re-running FoldX. Needs the HOLE progress file to already exist for the QC
column (same background, not the full-saturation one specifically -- QC
compares each label's CHAP radius against whichever HOLE run has it, and
41 populates results/tables/structural/pore_full_saturation_hole_progress.csv,
read here as well as the original background/literature HOLE outputs).
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
PDB_DIR = "data_processed/foldx_mutants/pore_full_saturation_pdbs"
WORK_ROOT = "data_processed/foldx_mutants/pore_full_saturation_chap"
ENUM_CSV = os.environ.get("ENUM_CSV", "results/tables/structural/pore_full_saturation_substitutions.csv")
PROGRESS_CSV = os.environ.get("PROGRESS_CSV", "results/tables/structural/pore_full_saturation_chap_progress.csv")
WT_CHAP_JSON = "data_processed/chap/kcna5_pore_repairedWT.json"
CHAP_IMAGE = "kcna5-chap:0.9.1"

INIT_POS_NM = np.array([13.2852, -13.28204, 13.31009])
DIR_VEC = np.array([0.01344207, 0.9998711, 0.008781797])
DIR_VEC = DIR_VEC / np.linalg.norm(DIR_VEC)
HYDRO_GRID = np.linspace(-3.8, 1.8, 300)


def load_hole_deltas():
    deltas = {}
    for path in [
        "results/tables/structural/pore_full_saturation_hole_progress.csv",
        "results/tables/structural/pore_background_hole_progress.csv",
    ]:
        if os.path.exists(path):
            with open(path, newline="") as f:
                for row in csv.DictReader(f):
                    if row["radius_delta_nm"] not in ("", "None"):
                        deltas.setdefault(row["label"], float(row["radius_delta_nm"]))
    return deltas


def ca_position(pdb_path, resno, chain="A"):
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
                if int(line[22:26]) == resno:
                    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return None


def already_done():
    if not os.path.exists(PROGRESS_CSV):
        return set()
    with open(PROGRESS_CSV, newline="") as f:
        return {row["label"] for row in csv.DictReader(f)}


_WT_HYDRO = None
_WT_RADIUS = None


def run_one(row):
    global _WT_HYDRO, _WT_RADIUS
    if _WT_HYDRO is None:
        wt = json.load(open(WT_CHAP_JSON))["pathwayProfile"]
        _WT_HYDRO = np.interp(HYDRO_GRID, wt["s"], wt["pfHydrophobicityMean"])
        _WT_RADIUS = (np.array(wt["s"]), np.array(wt["radiusMin"]))

    label = row["label"]
    position = int(row["position"])
    src_pdb = os.path.join(PDB_DIR, f"{label}.pdb")
    if not os.path.exists(src_pdb):
        return row, None, None

    work_dir = os.path.join(WORK_ROOT, label)
    os.makedirs(work_dir, exist_ok=True)
    hydro_max_delta, chap_radius_delta = None, None
    try:
        local_pdb = os.path.join(work_dir, "mutant.pdb")
        shutil.copy(src_pdb, local_pdb)

        abs_work = os.path.abspath(work_dir)
        chap_result = subprocess.run(
            ["docker", "run", "--rm", "--platform", "linux/amd64",
             "-v", f"{abs_work}:/out", "-w", "/out", CHAP_IMAGE,
             "-f", "/out/mutant.pdb", "-s", "/out/mutant.pdb",
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

        ca = ca_position(local_pdb, position)
        if ca is not None:
            s_axial = (ca / 10.0 - INIT_POS_NM) @ DIR_VEC
            m_s, m_r = np.array(prof["s"]), np.array(prof["radiusMin"])
            chap_radius_delta = float(np.interp(s_axial, m_s, m_r) - np.interp(s_axial, _WT_RADIUS[0], _WT_RADIUS[1]))
    except Exception:
        ## Same fix as 38/39/41: an uncaught docker-run timeout used to kill
        ## the whole worker pool instead of just failing this one row.
        hydro_max_delta, chap_radius_delta = None, None
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return row, hydro_max_delta, chap_radius_delta


def main():
    hole_deltas = load_hole_deltas()
    with open(ENUM_CSV, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["already_computed"] == "FALSE"]
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
        n_failed = 0
        with multiprocessing.Pool(N_WORKERS) as pool:
            for i, (row, hydro, chap_r) in enumerate(pool.imap_unordered(run_one, todo), 1):
                if hydro is None:
                    n_failed += 1  ## not written -- retried on next resume, see 38's rationale
                    continue
                hole_r = hole_deltas.get(row["label"])
                qc = abs(chap_r - hole_r) if (chap_r is not None and hole_r is not None) else None
                writer.writerow([row["label"], row["position"], row["wt_aa"], row["mut_aa"],
                                  hydro, chap_r, hole_r, qc])
                f.flush()
                if i % 20 == 0 or i == len(todo):
                    elapsed = time.time() - t_start
                    rate = elapsed / i
                    remaining = (len(todo) - i) * rate
                    print(f"[{i}/{len(todo)}] {row['label']} hydro_max_delta={hydro} qc={qc} "
                          f"({rate:.2f}s/mutation wall-clock, ~{remaining/60:.0f} min remaining)", flush=True)

    print(f"{n_failed} failed/timed-out this run (not recorded -- re-run to retry).", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
