#!/usr/bin/env python3
"""
Both interface orientations, side by side, for every observed (COSMIC)
T1-interface mutation -- not just whichever orientation happens to be more
extreme.

The "Chain" column in binding_ddg_expected_vs_observed.csv (inherited from
the Xu repo's AD/DA split) is NOT the PDB's own chain letters -- it's a tag
for which of the two FoldX AnalyseComplex runs (Binding_DDG_AD.csv /
Binding_DDG_DA.csv) a value came from. Those two runs correspond to the two
distinct local environments a T1-T1 interface residue sits in: KCNA5's
tetramer packs head-to-tail around a 4-fold ring, so a residue's contacts
looking "clockwise" to its neighbour are structurally different from its
contacts looking "counter-clockwise" to the other neighbour. Because the
mutation is present in the underlying gene, ALL FOUR subunit copies carry
it simultaneously in the real channel -- so both orientations are
perturbed at once in vivo, not just whichever one FoldX happened to score
higher. Quoting only the more extreme of the two (as informal discussion of
G128R had been doing: 8.5 kcal/mol from one orientation, 0.2 from the
other) overstates the case; this table reports both, explicitly, for every
mutation that has a value on both sides.
"""
import pandas as pd

RESULTS_DIR = "results"
ANALYSIS = "structural"

df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/binding_ddg_expected_vs_observed.csv")
observed = df[df["Group"] == "Observed"].copy()

wide = observed.pivot_table(index="Mutation", columns="Chain", values="DDG", aggfunc="first")
wide = wide.rename(columns={"A": "DDG_orientation_AD", "D": "DDG_orientation_DA"})
wide["both_orientations_tested"] = wide[["DDG_orientation_AD", "DDG_orientation_DA"]].notna().all(axis=1)
wide["max_abs_ddg"] = wide[["DDG_orientation_AD", "DDG_orientation_DA"]].abs().max(axis=1)
wide = wide.sort_values("max_abs_ddg", ascending=False).reset_index()

out_csv = f"{RESULTS_DIR}/tables/{ANALYSIS}/binding_ddg_both_orientations.csv"
wide.to_csv(out_csv, index=False)

print(f"n observed mutations: {len(wide)}")
print(f"n with both orientations tested: {wide['both_orientations_tested'].sum()}")
print()
print(wide.to_string(index=False))
print(f"\nSaved {out_csv}")
