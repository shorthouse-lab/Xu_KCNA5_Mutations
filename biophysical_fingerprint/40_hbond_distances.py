#!/usr/bin/env python3
"""
Two literature-implicated C-type-inactivation H-bond distances, measured
directly from each mutant tetramer PDB (kept by 38_pore_full_saturation_foldx.py
/ 39_regenerate_pdbs_for_hbonds.py in
data_processed/foldx_mutants/pore_full_saturation_pdbs/), for every mutation
across the full saturation set (1805 = 569 already-computed + 1236 new).

Biophysical rationale (full version in chat): C-type inactivation is
understood as the filter losing the rigid support that holds it in its
conductive geometry, not as something happening inside the filter itself.

  1. D485(OD1/OD2)-W472(NE1), intra-subunit (chain A only): the turret Asp
     "stapled" to the pore-helix-capping Trp immediately behind the filter.
     KCNA5 equivalent of Shaker's D447-W434 (Ple13) -- breaking it in Shaker
     (D447W) accelerated inactivation ~100-fold-class. Only defined when
     neither 472 nor 485 is itself mutated to something lacking the
     relevant atom.
  2. Y483(OH)-T477(OG1), inter-subunit ring (4-fold around the tetramer):
     Y483 is the filter's own Tyr reaching across to the NEXT chain's
     Thr477, coupling all 4 subunits' filter conformations together.
     Chain A (the only chain FoldX mutates here) plays both roles in this
     ring -- Y483 donor (paired with chain D's T477) and T477 acceptor
     (paired with chain B's Y483) -- so both instances are measured
     separately; a mutation at position 483 only breaks the first, a
     mutation at 477 only breaks the second, everything else affects
     neither atom directly (though could still shift the distance via
     local repacking).

WT reference values (KCNA5.pdb, verified in-session): D485-W472 = 2.87A,
Y483(A)-T477(D) = 2.90A, Y483(B)-T477(A) = 2.91A -- all squarely in H-bond
range, confirming these are real WT contacts, not spurious distances pulled
from the paper's text alone.
"""
import csv
import glob
import os

import numpy as np

PDB_DIR = "data_processed/foldx_mutants/pore_full_saturation_pdbs"
WT_PDB = "KCNA5.pdb"
OUT_CSV = "results/tables/structural/pore_full_saturation_hbond_distances.csv"


def atom_pos(pdb_path, chain, resno, atomname):
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == atomname:
                try:
                    if int(line[22:26]) == resno:
                        return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    continue
    return None


def hbond_features(pdb_path):
    d485_od1 = atom_pos(pdb_path, "A", 485, "OD1")
    d485_od2 = atom_pos(pdb_path, "A", 485, "OD2")
    w472_ne1 = atom_pos(pdb_path, "A", 472, "NE1")
    d485_w472 = None
    if w472_ne1 is not None:
        cands = [np.linalg.norm(o - w472_ne1) for o in (d485_od1, d485_od2) if o is not None]
        if cands:
            d485_w472 = float(min(cands))

    y483_a_oh = atom_pos(pdb_path, "A", 483, "OH")
    t477_d_og1 = atom_pos(pdb_path, "D", 477, "OG1")
    y483A_t477D = float(np.linalg.norm(y483_a_oh - t477_d_og1)) if (y483_a_oh is not None and t477_d_og1 is not None) else None

    y483_b_oh = atom_pos(pdb_path, "B", 483, "OH")
    t477_a_og1 = atom_pos(pdb_path, "A", 477, "OG1")
    y483B_t477A = float(np.linalg.norm(y483_b_oh - t477_a_og1)) if (y483_b_oh is not None and t477_a_og1 is not None) else None

    return d485_w472, y483A_t477D, y483B_t477A


def main():
    wt_d485_w472, wt_y483A_t477D, wt_y483B_t477A = hbond_features(WT_PDB)
    print(f"WT: D485-W472={wt_d485_w472:.3f} Y483(A)-T477(D)={wt_y483A_t477D:.3f} Y483(B)-T477(A)={wt_y483B_t477A:.3f}")

    pdbs = sorted(glob.glob(os.path.join(PDB_DIR, "*.pdb")))
    print(f"{len(pdbs)} mutant PDBs found")

    rows = []
    n_missing_atom = 0
    for pdb_path in pdbs:
        label = os.path.splitext(os.path.basename(pdb_path))[0]
        d, yA, yB = hbond_features(pdb_path)
        if d is None or yA is None or yB is None:
            n_missing_atom += 1
        rows.append(dict(
            label=label,
            d485_w472_dist=d,
            y483A_t477D_dist=yA,
            y483B_t477A_dist=yB,
            d485_w472_delta=(d - wt_d485_w472) if d is not None else None,
            y483A_t477D_delta=(yA - wt_y483A_t477D) if yA is not None else None,
            y483B_t477A_delta=(yB - wt_y483B_t477A) if yB is not None else None,
        ))

    print(f"{n_missing_atom} labels missing at least one atom (mutation removed a bond-partner side chain, e.g. at 472/477/483/485 itself)")

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {OUT_CSV}")


if __name__ == "__main__":
    main()
