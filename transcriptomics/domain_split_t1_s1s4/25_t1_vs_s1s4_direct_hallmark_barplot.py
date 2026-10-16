#!/usr/bin/env python3
"""
Direct T1-vs-S1S4 contrast, Hallmark only -- same plot style as
24_allmutant_pathway_barplot.py's make_single_barplot (matplotlib, orange/
green by NES sign, seaborn despine) so this figure reads consistently with
the WT-vs-AllMutant Hallmark barchart in the first main figure. Source data
is run_gsea_t1_vs_s1s4_direct_combined.R's already-filtered
fgsea_T1_vs_S1S4_combined_solid_sig.csv (FDR q<0.25 AND leading edge >=5
genes -- the direct contrast is underpowered (n=11 vs n=11) relative to the
vs-WT contrasts, hence the GSEA exploratory 0.25 threshold rather than 0.05;
see that script's header for the leading-edge thin-evidence rationale),
restricted here to collection=="hallmark".
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RESULTS_DIR = "results"
ANALYSIS = "analysis3_domain_split"


def clean_pathway_label(name):
    label = name.replace("HALLMARK_", "", 1).replace("_", " ")
    return label if len(label) <= 50 else label[:47] + "..."


def make_single_barplot(df, title, out_name, xlabel):
    top = df.sort_values("NES").reset_index(drop=True)
    n = len(top)

    fig, ax = plt.subplots(figsize=(9, max(6.5, 0.45 * n + 2)))
    # Matches the venn_genes/pathways_*_T1_vs_S1S4 palette: NES<0 = higher in
    # S1S4 (purple), NES>0 = higher in T1 (orange).
    colors = ["#7570b3" if nes < 0 else "#d95f02" for nes in top["NES"]]
    bars = ax.barh(top["pathway"], top["NES"], color=colors)

    xmax = top["NES"].abs().max()
    ax.set_xlim(-xmax * 1.35, xmax * 1.35)
    for bar, padj, n_le in zip(bars, top["padj"], top["n_leading_edge"]):
        w = bar.get_width()
        label = f"FDR q={padj:.1e}" if padj < 0.001 else f"FDR q={padj:.3f}"
        label += f"  ({n_le} genes)"
        ax.text(-(1 if w == 0 else (1 if w > 0 else -1)) * xmax * 0.03, bar.get_y() + bar.get_height() / 2, label,
                va="center", ha="right" if w >= 0 else "left", fontsize=12, color="dimgrey")

    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=13.5)
    ax.set_title(title, fontsize=16, fontweight="bold")
    sns.despine(ax=ax)
    fig.tight_layout()

    out_png = f"{RESULTS_DIR}/figures/{ANALYSIS}/{out_name}.png"
    out_pdf = f"{RESULTS_DIR}/figures/{ANALYSIS}/{out_name}.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"Saved {out_png} and {out_pdf}")


df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/fgsea_T1_vs_S1S4_combined_solid_sig.csv")
df = df[df["collection"] == "hallmark"].copy()
print(f"Hallmark pathways directly distinguishing T1 from S1S4 (FDR q<0.25, leading edge >=5): {len(df)}")
df["pathway"] = df["pathway"].apply(clean_pathway_label)

make_single_barplot(
    df,
    "Hallmark pathways directly distinguishing T1 from S1S4 (FDR q<0.25)",
    "fgsea_T1_vs_S1S4_hallmark_barplot",
    "NES (direct T1 vs S1S4 GSEA, Hallmark)",
)
