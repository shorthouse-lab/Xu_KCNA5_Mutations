#!/usr/bin/env python3
"""
Single-panel version of 22_hypoxia_domain_split_seaborn.py's 3-method figure
-- same box+stripplot style, same domain groups/palette/pairwise-bracket
logic, but just ONE hypoxia-scoring method instead of all three side by
side. Method choice: hif_panel_score (the literature HIF-target gene panel,
z-scored), because it's the one with the strongest T1-vs-S1S4 signal of the
three (nominal p=0.025, vs 0.076 for PROGENy and 0.131 for Hallmark ssGSEA
-- see hypoxia_domain_pairwise_all_models.csv), and HALLMARK_HYPOXIA itself
also showed up as a borderline pathway in the direct T1-vs-S1S4 GSEA
(FDR q=0.154), so this cross-checks that finding with an independent,
non-GSEA scoring method on the same biology.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RESULTS_DIR = "results"
ANALYSIS = "analysis4_network_immune"
P_ANNOTATE_THRESHOLD = 0.1
METHOD_COL = "hif_panel_score"
METHOD_LABEL = "HIF-target panel (z-score)"

GROUP_ORDER = ["T1_mutant", "S1S4_mutant", "WT"]  # Other_domain_mutant and S5S6_mutant (n=1) dropped -- neither T1 nor S1S4
PALETTE = {"T1_mutant": "#d95f02", "S1S4_mutant": "#7570b3", "WT": "#D9D9D9"}  # matches the venn_genes/pathways_*_T1_vs_S1S4 palette
LABELS = {"T1_mutant": "T1", "S1S4_mutant": "S1-S4", "WT": "WT"}
X_POS = {g: i for i, g in enumerate(GROUP_ORDER)}

df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/hypoxia_scores_by_domain_group.csv")
df = df[df["domain_group"].isin(GROUP_ORDER)].copy()

pairwise = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/hypoxia_domain_pairwise_all_models.csv")

df["domain_group"] = pd.Categorical(df["domain_group"], categories=GROUP_ORDER, ordered=True)


def p_to_text(p):
    return f"p={p:.3f}"


def add_pairwise_brackets(ax, col, data_max):
    # Raw (uncorrected) nominal p-value, both for selecting which pairs to
    # annotate and for the label itself -- this is the one method/figure
    # being shown, not a multi-method/multi-group family to correct across.
    sig_pairs = pairwise[(pairwise["score"] == col) & (pairwise["p_value"] < P_ANNOTATE_THRESHOLD) &
                          pairwise["group1"].isin(GROUP_ORDER) & pairwise["group2"].isin(GROUP_ORDER)].copy()
    if sig_pairs.empty:
        return
    sig_pairs["span"] = (sig_pairs["group1"].map(X_POS) - sig_pairs["group2"].map(X_POS)).abs()
    sig_pairs = sig_pairs.sort_values("span")

    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    step = y_range * 0.09
    y = data_max + y_range * 0.06
    for _, row in sig_pairs.iterrows():
        x1, x2 = sorted([X_POS[row["group1"]], X_POS[row["group2"]]])
        ax.plot([x1, x1, x2, x2], [y, y + step * 0.15, y + step * 0.15, y], lw=1.2, color="black")
        ax.text((x1 + x2) / 2, y + step * 0.2, p_to_text(row["p_value"]),
                ha="center", va="bottom", fontsize=9.5)
        y += step
    ax.set_ylim(ax.get_ylim()[0], y + step * 0.5)


fig, ax = plt.subplots(figsize=(7.5, 6))

sns.boxplot(data=df, x="domain_group", y=METHOD_COL, order=GROUP_ORDER, hue="domain_group",
            palette=PALETTE, showfliers=False, width=0.55, linewidth=1.5, legend=False, ax=ax)
sns.stripplot(data=df, x="domain_group", y=METHOD_COL, order=GROUP_ORDER, color="black",
              alpha=0.3, size=3.5, jitter=0.15, ax=ax)
ax.set_xlabel("")
ax.set_ylabel(METHOD_LABEL, fontsize=14)
ax.set_xticks(range(len(GROUP_ORDER)))
ax.set_xticklabels([LABELS[g] for g in GROUP_ORDER], fontsize=13)
ax.tick_params(axis="y", labelsize=12)
add_pairwise_brackets(ax, METHOD_COL, df[METHOD_COL].max())
sns.despine(ax=ax)

n_counts = df.groupby("domain_group", observed=True).size()
fig.suptitle("Hypoxia score by KCNA5 structural-domain group", fontsize=15, y=0.97)
fig.text(0.5, 0.905, " / ".join(f"{LABELS[g]} n={n_counts.get(g, 0)}" for g in GROUP_ORDER),
         ha="center", fontsize=11, color="dimgrey")
fig.text(0.5, 0.01,
         "p-value: MSI+TMB-adjusted linear model, raw (uncorrected)",
         ha="center", fontsize=9, color="dimgrey", style="italic")

fig.tight_layout(rect=(0, 0.08, 1, 0.86))

out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/hypoxia_domain_split_hif_panel_single.png"
out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/hypoxia_domain_split_hif_panel_single.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved {out_png} and {out_pdf}")
