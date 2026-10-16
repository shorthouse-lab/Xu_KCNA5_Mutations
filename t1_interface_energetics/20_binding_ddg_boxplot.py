#!/usr/bin/env python3
"""
Re-plot the Xu repo's binding DDG analysis (binding_ddg/ddg_analysis.ipynb)
as two boxplots instead of its overlaid histogram/KDE:
  - "Expected": DDG for every possible substitution FoldX tested at every T1
    inter-subunit interface residue (Binding_DDG_AD.csv + _DA.csv), i.e. the
    background distribution if mutations hit the interface with no selection.
  - "Observed": the subset of those substitutions that are actual COSMIC COAD
    mutations (KCNA5_T1_mutations.csv, the Xu repo's own binding_ddg input
    list -- includes R121H/R121C, just outside the T1_DOMAIN_START/END
    (128-216) boundary used elsewhere in this project, since they're still
    on the FoldX-tested interface) AND sit on an interface residue.

Mutation-code reformatting (FoldX "FE196D" -> "F196D") is reused verbatim
from ddg_analysis.ipynb rather than re-derived.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

RESULTS_DIR = "results"
ANALYSIS = "structural"
XU_DIR = "Xu_KCNA5_Mutations-main/binding_ddg"

# =============================================================================
# A. Load + reformat, reused from ddg_analysis.ipynb
# =============================================================================
def load_ddg(path, chain_label):
    df = pd.read_csv(path, sep=",")
    df["Mutation"] = [m[0] + str(int(m[2:-1])) + m[-1] for m in df["Mutation"]]
    df["Chain"] = chain_label
    return df


ddg_df_AD = load_ddg(f"{XU_DIR}/results/Binding_DDG_AD.csv", "A")
ddg_df_DA = load_ddg(f"{XU_DIR}/results/Binding_DDG_DA.csv", "D")
ddg_df = pd.concat([ddg_df_AD, ddg_df_DA], ignore_index=True)
print(f"Loaded {len(ddg_df)} interface substitution DDG values "
      f"(chain A: {len(ddg_df_AD)}, chain D: {len(ddg_df_DA)})")

cosmic_df = pd.read_csv(f"{XU_DIR}/data/KCNA5_T1_mutations.csv")
n_t1_unique = cosmic_df["Mutation"].nunique()
print(f"Loaded {n_t1_unique} unique COSMIC COAD T1-domain mutations")

ddg_df_cosmic = ddg_df[ddg_df["Mutation"].isin(cosmic_df["Mutation"])]
print(f"Observed (on-interface COSMIC mutations) found in interface DDG data: "
      f"{ddg_df_cosmic['Mutation'].nunique()} of {n_t1_unique} COSMIC T1-domain mutations "
      f"({len(ddg_df_cosmic)} rows across both chain orientations)")

# =============================================================================
# B. Assemble plotting frame + Mann-Whitney test
# =============================================================================
plot_df = pd.concat([
    ddg_df.assign(Group="Expected"),
    ddg_df_cosmic.assign(Group="Observed"),
], ignore_index=True)
plot_df["Group"] = pd.Categorical(plot_df["Group"], categories=["Expected", "Observed"])

stat, p = mannwhitneyu(ddg_df["DDG"], ddg_df_cosmic["DDG"], alternative="two-sided")
print(f"Expected vs observed DDG, Mann-Whitney U={stat:.1f}, p={p:.3g} "
      f"(median expected={ddg_df['DDG'].median():.3f}, median observed={ddg_df_cosmic['DDG'].median():.3f})")

out_dir_tables = f"{RESULTS_DIR}/tables/{ANALYSIS}"
plot_df[["Mutation", "Chain", "WT_energy", "MUT_energy", "DDG", "Group"]].to_csv(
    f"{out_dir_tables}/binding_ddg_expected_vs_observed.csv", index=False
)
pd.DataFrame({"U_statistic": [stat], "p_value": [p]}).to_csv(
    f"{out_dir_tables}/binding_ddg_mannwhitney_test.csv", index=False
)

# =============================================================================
# C. Two boxplots: Expected vs Observed
# =============================================================================
PALETTE = {"Expected": "#2a78d6", "Observed": "#eb6834"}

fig, ax = plt.subplots(figsize=(7, 6.5))
sns.boxplot(data=plot_df, x="Group", y="DDG", hue="Group", palette=PALETTE,
            showfliers=False, width=0.5, linewidth=1.5, legend=False, ax=ax)
sns.stripplot(data=plot_df, x="Group", y="DDG", color="black", alpha=0.3, size=3,
              jitter=0.15, ax=ax)

ax.set_xlabel("")
ax.set_ylabel("Binding ΔΔG (kcal/mol)", fontsize=14)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"Expected\n(n={len(ddg_df)})", f"Observed\n(n={len(ddg_df_cosmic)})"])
fig.suptitle("Interface binding ΔΔG: all possible substitutions vs.\nCOSMIC-observed mutations",
             fontsize=13, y=0.99)
ax.set_title(f"Mann-Whitney p = {p:.3g}", fontsize=10)
ax.tick_params(axis="both", labelsize=12)
sns.despine(ax=ax)
fig.tight_layout(rect=(0, 0, 1, 0.90))

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/binding_ddg_expected_vs_observed_boxplot.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/binding_ddg_expected_vs_observed_boxplot.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
