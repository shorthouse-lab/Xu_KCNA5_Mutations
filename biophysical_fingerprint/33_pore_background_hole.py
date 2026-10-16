#!/usr/bin/env python3
"""
HOLE2 pore-radius change (radius_delta_nm) for the 549 background
substitutions, to sit alongside the FoldX energy decomposition (30/31/32)
in the PCA/permutation analysis -- adding back one of the two dimensions
originally scoped out (see chat) now that this has been piloted and
validated end-to-end.

Pipeline per mutation, all disposable (no per-mutation FoldX/HOLE output
kept beyond this run, same rationale as 30_pore_background_foldx.py):
  1. FoldX BuildModel (same convention as 30_pore_background_foldx.py,
     but this time the mutant tetramer PDB is kept, not deleted, since
     HOLE needs it) -- re-run rather than reusing 30's already-computed
     ddG, since 30 didn't keep the structures.
  2. Trim to the pore-forming module (chains A-D, residues 425-527) --
     replicates the boundary used for the original 13 mutations
     (data_processed/hole/mutants/*_pore_module.pdb, verified in-session:
     3194 atoms for a known example, chain A only carries the mutation,
     matching FoldX's heteromeric convention). Plain fixed-column PDB
     line filtering -- no bio3d/biopython dependency needed for this.
  3. HOLE2 via Docker (kcna5-hole:latest, run via stdin -- confirmed in
     pilot that passing the .inp as a CLI arg does NOT work, HOLE reads
     control input from stdin only), same CPOINT/CVECT as the validated
     WT axis (13_chap_pore_analysis.R), same SAMPLE/ENDRAD/MCDISP/MCSTEP
     as the original 13's config (data_processed/hole/mutants/*.inp).
  4. radius_delta_nm: interpolate the mutant's own HOLE centerline at the
     mutated residue's own axial position (CA projected onto the
     validated pore axis), same as WT's centerline at that position, take
     the difference -- identical method to
     21_pore_radius_hydrophobicity_mutants.py's pore_radius_nm_centerline
     column, verified in pilot against the two already-known observed
     values (S449P ~0, G504R -0.008) before running this at scale.

Pilot timing (5 mutations, serial): FoldX ~13s/mutation (unparallelised
single run), HOLE ~1.5s/mutation including Docker container startup.
HOLE is not the bottleneck here -- FoldX is, same as before, so this
reuses the same N_WORKERS process-pool pattern as
30_pore_background_foldx.py.
"""
import csv
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
WORK_ROOT = "data_processed/foldx_mutants/pore_background_hole"
ENUM_CSV = os.environ.get("ENUM_CSV", "results/tables/structural/pore_background_substitutions.csv")
PROGRESS_CSV = os.environ.get("PROGRESS_CSV", "results/tables/structural/pore_background_hole_progress.csv")
WT_SPH = "data_processed/hole/mutants/wt_hole.sph"

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


def load_enum_rows():
    with open(ENUM_CSV, newline="") as f:
        return list(csv.DictReader(f))


def already_done():
    if not os.path.exists(PROGRESS_CSV):
        return set()
    with open(PROGRESS_CSV, newline="") as f:
        return {row["label"] for row in csv.DictReader(f)}


_WT_S, _WT_R = None, None  # loaded once per worker process, not per mutation


def run_one(row):
    global _WT_S, _WT_R
    if _WT_S is None:
        _WT_S, _WT_R = centerline(WT_SPH)

    label = row["label"]
    wt_aa, position, mut_aa = row["wt_aa"], int(row["position"]), row["mut_aa"]
    work_dir = os.path.join(WORK_ROOT, label)
    os.makedirs(work_dir, exist_ok=True)

    radius_delta_nm = None
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
            return row, None

        pore_pdb = os.path.join(work_dir, "pore_module.pdb")
        trim_pore_module(mutant_pdb, pore_pdb)

        sph_path = os.path.join(work_dir, "hole.sph")
        inp = INP_TEMPLATE.format(pdb=f"/out/{label}/pore_module.pdb", sph=f"/out/{label}/hole.sph")
        hole_result = subprocess.run(
            ["docker", "run", "--rm", "-i", "-v", f"{os.path.abspath(WORK_ROOT)}:/out",
             "kcna5-hole:latest", "hole"],
            input=inp, capture_output=True, text=True, timeout=60,
        )
        if not os.path.exists(sph_path):
            sys.stderr.write(f"[FAILED HOLE] {label}: {hole_result.stdout[-300:]}\n")
            return row, None

        ca = ca_position(pore_pdb, position)
        if ca is None:
            sys.stderr.write(f"[FAILED] {label}: CA not found in trimmed pore module\n")
            return row, None
        s_axial = (ca / 10.0 - INIT_POS_NM) @ DIR_VEC
        m_s, m_r = centerline(sph_path)
        mutant_radius = np.interp(s_axial, m_s, m_r)
        wt_radius = np.interp(s_axial, _WT_S, _WT_R)
        radius_delta_nm = float(mutant_radius - wt_radius)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return row, radius_delta_nm


def main():
    rows = load_enum_rows()
    done = already_done()
    todo = [r for r in rows if r["label"] not in done]
    print(f"{len(rows)} total, {len(done)} already done, {len(todo)} to run ({N_WORKERS} workers)", flush=True)

    write_header = not os.path.exists(PROGRESS_CSV)
    with open(PROGRESS_CSV, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["label", "position", "wt_aa", "mut_aa", "radius_delta_nm"])
            f.flush()

        t_start = time.time()
        with multiprocessing.Pool(N_WORKERS) as pool:
            for i, (row, delta) in enumerate(pool.imap_unordered(run_one, todo), 1):
                writer.writerow([row["label"], row["position"], row["wt_aa"], row["mut_aa"], delta])
                f.flush()
                if i % 10 == 0 or i == len(todo):
                    elapsed = time.time() - t_start
                    rate = elapsed / i
                    remaining = (len(todo) - i) * rate
                    print(f"[{i}/{len(todo)}] {row['label']} radius_delta_nm={delta} "
                          f"({rate:.2f}s/mutation wall-clock, ~{remaining/60:.0f} min remaining)", flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
