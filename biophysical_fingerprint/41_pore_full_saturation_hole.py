#!/usr/bin/env python3
"""
HOLE2 pore-radius change for the full saturation set (1805 = 19 non-WT AAs x
95 positions, 433-527) -- extends 33_pore_background_hole.py's already-
validated method (trim to pore module 425-527/chains A-D, HOLE via Docker
stdin, radius at the mutated residue's own axial position vs. WT centerline)
from the 549 codon-reachable background to every possible substitution.

Unlike 33, this reuses the mutant tetramer PDBs already on disk
(data_processed/foldx_mutants/pore_full_saturation_pdbs/<label>.pdb, saved
by 38_pore_full_saturation_foldx.py / 39_regenerate_pdbs_for_hbonds.py)
rather than re-running FoldX -- BuildModel is deterministic for a given
substitution (verified repeatedly this session: identical ddG across
independent re-runs), so there's no reason to redo it a third time. This
also means only the ~1236 rows that don't already have a HOLE value (549
background + 9 literature already computed in 33/the ad hoc literature
script) need running here.
"""
import csv
import multiprocessing
import os

import numpy as np

N_WORKERS = 8
PDB_DIR = "data_processed/foldx_mutants/pore_full_saturation_pdbs"
WORK_ROOT = "data_processed/foldx_mutants/pore_full_saturation_hole"
ENUM_CSV = os.environ.get("ENUM_CSV", "results/tables/structural/pore_full_saturation_substitutions.csv")
PROGRESS_CSV = os.environ.get("PROGRESS_CSV", "results/tables/structural/pore_full_saturation_hole_progress.csv")
WT_HOLE_SPH = "data_processed/hole/mutants/wt_hole.sph"

INIT_POS_NM = np.array([13.2852, -13.28204, 13.31009])
DIR_VEC = np.array([0.01344207, 0.9998711, 0.008781797])
DIR_VEC = DIR_VEC / np.linalg.norm(DIR_VEC)

INP_TEMPLATE = """COORD {pdb}
RADIUS /opt/conda/share/hole2/rad/simple.rad
SPHPDB {sph}
CPOINT 132.852 -132.8204 133.1009
CVECT 0.01344207 0.9998711 0.008781797
SAMPLE 0.2
ENDRAD 10.0
MCDISP 0.1
MCSTEP 500
"""


def load_sph(path):
    xyz, radius = [], []
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            radius.append(float(line[54:60]))
    return np.array(xyz), np.array(radius)


def centerline(sph_path):
    xyz, radius = load_sph(sph_path)
    xyz_nm = xyz / 10.0
    s_axial = (xyz_nm - INIT_POS_NM) @ DIR_VEC
    order = np.argsort(s_axial)
    return s_axial[order], (radius / 10.0)[order]


def ca_position(pdb_path, resno, chain="A"):
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
                if int(line[22:26]) == resno:
                    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return None


def trim_pore_module(src_pdb, dst_pdb, chains=("A", "B", "C", "D"), resno_range=(425, 527)):
    lo, hi = resno_range
    with open(src_pdb) as fin, open(dst_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                try:
                    resno = int(line[22:26])
                except ValueError:
                    continue
                if chain in chains and lo <= resno <= hi:
                    fout.write(line)
            elif line.startswith(("TER", "END")):
                fout.write(line)


def already_done():
    if not os.path.exists(PROGRESS_CSV):
        return set()
    with open(PROGRESS_CSV, newline="") as f:
        return {row["label"] for row in csv.DictReader(f)}


_WT_S, _WT_R = None, None


def run_one(row):
    global _WT_S, _WT_R
    if _WT_S is None:
        _WT_S, _WT_R = centerline(WT_HOLE_SPH)

    label = row["label"]
    position = int(row["position"])
    src_pdb = os.path.join(PDB_DIR, f"{label}.pdb")
    if not os.path.exists(src_pdb):
        return row, None  # PDB not yet generated (shouldn't happen if 38/39 ran to completion)

    work_dir = os.path.join(WORK_ROOT, label)
    os.makedirs(work_dir, exist_ok=True)
    radius_delta_nm = None
    try:
        pore_pdb = os.path.join(work_dir, "pore_module.pdb")
        trim_pore_module(src_pdb, pore_pdb)

        sph_path = os.path.join(work_dir, "hole.sph")
        inp = INP_TEMPLATE.format(pdb=f"/out/{label}/pore_module.pdb", sph=f"/out/{label}/hole.sph")
        import subprocess
        subprocess.run(
            ["docker", "run", "--rm", "-i", "-v", f"{os.path.abspath(WORK_ROOT)}:/out",
             "kcna5-hole:latest", "hole"],
            input=inp, capture_output=True, text=True, timeout=60,
        )
        if not os.path.exists(sph_path):
            return row, None

        ca = ca_position(pore_pdb, position)
        if ca is None:
            return row, None
        s_axial = (ca / 10.0 - INIT_POS_NM) @ DIR_VEC
        m_s, m_r = centerline(sph_path)
        radius_delta_nm = float(np.interp(s_axial, m_s, m_r) - np.interp(s_axial, _WT_S, _WT_R))
    except Exception:
        ## Same fix as 38/39: an uncaught docker-run timeout (or any other
        ## exception) here used to kill the whole worker pool instead of
        ## just failing this one row.
        radius_delta_nm = None
    finally:
        import shutil
        shutil.rmtree(work_dir, ignore_errors=True)

    return row, radius_delta_nm


def main():
    with open(ENUM_CSV, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["already_computed"] == "FALSE"]
    done = already_done()
    todo = [r for r in rows if r["label"] not in done]
    print(f"{len(rows)} total, {len(done)} already done, {len(todo)} to run ({N_WORKERS} workers)", flush=True)

    write_header = not os.path.exists(PROGRESS_CSV)
    with open(PROGRESS_CSV, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["label", "position", "wt_aa", "mut_aa", "radius_delta_nm"])
            f.flush()

        import time
        t_start = time.time()
        n_failed = 0
        with multiprocessing.Pool(N_WORKERS) as pool:
            for i, (row, delta) in enumerate(pool.imap_unordered(run_one, todo), 1):
                if delta is None:
                    n_failed += 1  ## not written -- retried on next resume, see 38's rationale
                    continue
                writer.writerow([row["label"], row["position"], row["wt_aa"], row["mut_aa"], delta])
                f.flush()
                if i % 20 == 0 or i == len(todo):
                    elapsed = time.time() - t_start
                    rate = elapsed / i
                    remaining = (len(todo) - i) * rate
                    print(f"[{i}/{len(todo)}] {row['label']} radius_delta_nm={delta} "
                          f"({rate:.2f}s/mutation wall-clock, ~{remaining/60:.0f} min remaining)", flush=True)

    print(f"{n_failed} failed/timed-out this run (not recorded -- re-run to retry).", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
