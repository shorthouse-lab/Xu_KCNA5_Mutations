#!/usr/bin/env python3
"""
Single-contrast (AllMutant vs WT) pathway bar charts -- the Part 1 pathway
figures, companion to volcano_AllMutant_vs_WT.png. Unlike the T1-vs-S1S4
paired bar charts (23_top10_paired_barplots.py), there's only one contrast
here, so each is a plain horizontal NES bar chart: ALL pathways passing
FDR q<0.05 (padj IS the BH-adjusted FDR q-value, not a different statistic),
signed NES so direction reads directly off the bar. Built for all 4 GSEA
collections (Hallmark, Reactome, KEGG, Oncogenic), same as the T1-vs-S1S4
pathway comparison in Part 2. Large collections are capped at max_show for
legibility (full set always saved to a companion *_all_passing.csv).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RESULTS_DIR = "results"
ANALYSIS = "analysis3_transcriptomic"
FDR_THRESHOLD = 0.05
MAX_SHOW = 50


def clean_pathway_label(name, prefix):
    label = name.replace(f"{prefix}_", "", 1).replace("_", " ") if prefix else name.replace("_", " ")
    return label if len(label) <= 50 else label[:47] + "..."


def make_single_barplot(df, title, out_name, xlabel, fdr_threshold=FDR_THRESHOLD, max_show=MAX_SHOW):
    n_before = len(df)
    df = df[df["padj"] < fdr_threshold].copy()
    n_passing = len(df)
    print(f"  [{out_name}] FDR q<{fdr_threshold} filter: kept {n_passing}/{n_before} pathways")
    if n_passing == 0:
        print(f"  [{out_name}] no pathways pass the FDR threshold, skipping figure")
        return

    table_path = f"{RESULTS_DIR}/tables/{ANALYSIS}/{out_name}_all_passing.csv"
    df.sort_values("padj").to_csv(table_path, index=False)

    if n_passing > max_show:
        print(f"  [{out_name}] capping display at top {max_show} of {n_passing} by FDR q (all {n_passing} saved to {table_path})")
        df = df.sort_values("padj").head(max_show)
        title = f"{title}\n(top {max_show} of {n_passing} pathways at FDR q<{fdr_threshold}, by FDR q)"
    top = df.sort_values("NES").reset_index(drop=True)
    n = len(top)

    # Per-bar FDR-value text collides with the y-axis labels once there are
    # many rows -- see 23_top10_paired_barplots.py for the same fix. Above
    # ~20 rows, drop the inline annotation; exact values are in the CSV.
    show_annotations = n <= 20
    if not show_annotations:
        title = f"{title}\n(exact FDR q-values for all shown pathways in {out_name}_all_passing.csv)"

    fig, ax = plt.subplots(figsize=(9, max(6.5, 0.45 * n + 2)))
    colors = ["#d95f02" if nes < 0 else "#1b9e77" for nes in top["NES"]]
    bars = ax.barh(top["pathway"], top["NES"], color=colors)

    xmax = top["NES"].abs().max()
    ax.set_xlim(-xmax * 1.35, xmax * 1.35)
    if show_annotations:
        for bar, padj in zip(bars, top["padj"]):
            w = bar.get_width()
            label = f"FDR q={padj:.1e}" if padj < 0.001 else f"FDR q={padj:.3f}"
            # Labels sit just past the zero line, on the side opposite the bar's
            # own extent -- the far end of a bar can run into the y-axis tick
            # labels when most bars share one sign, but the area beyond zero is
            # always empty regardless of direction mix.
            ax.text(-np.sign(w if w != 0 else 1) * xmax * 0.03, bar.get_y() + bar.get_height() / 2, label,
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


COLLECTIONS = [
    ("hallmark", "HALLMARK", "Hallmark pathways at FDR q<0.05: All KCNA5 mutants vs WT", "NES (fgsea, Hallmark)"),
    ("reactome", "REACTOME", "Reactome pathways at FDR q<0.05: All KCNA5 mutants vs WT", "NES (fgsea, Reactome)"),
    ("kegg", "KEGG", "KEGG pathways at FDR q<0.05: All KCNA5 mutants vs WT", "NES (fgsea, KEGG)"),
    ("oncogenic", None, "Oncogenic Signatures (C6) at FDR q<0.05: All KCNA5 mutants vs WT", "NES (fgsea, Oncogenic)"),
]

for coll, prefix, title, xlabel in COLLECTIONS:
    df = pd.read_csv(f"{RESULTS_DIR}/tables/{ANALYSIS}/fgsea_AllMutant_vs_WT_{coll}.csv")
    df = df[~df["padj"].isna()].copy()
    df["pathway"] = df["pathway"].apply(lambda n: clean_pathway_label(n, prefix))
    make_single_barplot(df, title, f"{coll}_allmutant_vs_wt_barplot", xlabel)
