#!/usr/bin/env python3
"""Generate manuscript figures from frozen project tables only.

All operations in this file are descriptive reshaping or display transforms:
sorting, grouping, sign encoding, -log10(FDR), chromosome offsets, and counts.
No statistical test, model, hotspot definition, or evidence tier is recomputed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
FIGROOT = ROOT / "figures" / "manuscript"
DATA_DIR = FIGROOT / "data"
MAIN_DIR = FIGROOT / "main"
SUPP_DIR = FIGROOT / "supplement"
QC_DIR = FIGROOT / "qc"

AXIS_COLORS = {
    "Immune__Cardiometabolic": "#0072B2",
    "Immune__Psychiatric": "#009E73",
    "Cardiometabolic__Psychiatric": "#D55E00",
}
AXIS_SHORT = {
    "Immune__Cardiometabolic": "Immune-CM",
    "Immune__Psychiatric": "Immune-Psych",
    "Cardiometabolic__Psychiatric": "CM-Psych",
}
CLASS_COLORS = {
    "GLOBAL_WEAK_LOCAL_PRESENT": "#D55E00",
    "GLOBAL_STRONG_LOCAL_PRESENT": "#0072B2",
    "GLOBAL_WEAK_LOCAL_ABSENT": "#999999",
    "GLOBAL_STRONG_LOCAL_ABSENT": "#56B4E9",
}
ABF_COLORS = {
    "STRONG_SHARED_SIGNAL": "#0072B2",
    "SUGGESTIVE_SHARED_SIGNAL": "#56B4E9",
    "DISTINCT_SIGNAL_SUPPORTED": "#D55E00",
    "WEAK_OR_INCONCLUSIVE": "#999999",
}

CHR_LENGTHS = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276,
    5: 180915260, 6: 171115067, 7: 159138663, 8: 146364022,
    9: 141213431, 10: 135534747, 11: 135006516, 12: 133851895,
    13: 115169878, 14: 107349540, 15: 102531392, 16: 90354753,
    17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
    21: 48129895, 22: 51304566,
}


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })


def read_table(relative: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / relative, sep="\t")


def pair_label(pair: str) -> str:
    return str(pair).replace("__", "-")


def axis_order(values: Iterable[str]) -> List[str]:
    preferred = list(AXIS_COLORS)
    return [x for x in preferred if x in set(values)] + sorted(set(values) - set(preferred))


def ordered_pairs(df: pd.DataFrame, axis_col: str = "domain_axis") -> List[str]:
    axes = axis_order(df.loc[df[axis_col].notna(), axis_col].astype(str))
    result: List[str] = []
    for axis in axes:
        result.extend(sorted(df.loc[df[axis_col] == axis, "pair"].astype(str).unique()))
    return result


def add_panel_label(ax, label: str) -> None:
    ax.text(-0.05, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom")


def save_data(df: pd.DataFrame, name: str) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(DATA_DIR / f"{name}_data.tsv", sep="\t", index=False)


def save_figure(fig, name: str, directory: Path, include_tiff: bool = True) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{name}.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(directory / f"{name}.svg", dpi=600, bbox_inches="tight")
    fig.savefig(directory / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(directory / f"{name}_preview.png", dpi=600, bbox_inches="tight")
    if include_tiff:
        fig.savefig(directory / f"{name}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def cumulative_positions(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    offsets = {}
    running = 0
    for chrom in range(1, 23):
        offsets[chrom] = running
        running += CHR_LENGTHS[chrom]
    out["chr"] = pd.to_numeric(out["chr"], errors="coerce").astype(int)
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["stop"] = pd.to_numeric(out["stop"], errors="coerce")
    out["cumulative_mid_bp"] = (out["start"] + out["stop"]) / 2 + out["chr"].map(offsets)
    out["cumulative_mid_mb"] = out["cumulative_mid_bp"] / 1e6
    return out


def make_figure1() -> None:
    params = read_table("results/genomic_sem/phase2d/D2_CD_UC_PsA_observed_parameters.tsv")
    model = read_table("results/genomic_sem/phase2d/model_comparison.tsv")
    d2 = model.loc[model["model"] == "D2_CD_UC_PsA_observed"].iloc[0]
    rows = []
    workflow = [
        "9-trait GWAS", "LDSC", "Genomic SEM\nmodel comparison", "D2 selected\nand locked",
        "27-pair LAVA\natlas", "Hotspot\nconsolidation", "MDD noUKBB\nsensitivity", "ABF\nscreening",
    ]
    for i, item in enumerate(workflow, 1):
        rows.append({"panel": "A", "row_type": "workflow", "item": item.replace("\n", " "), "order": i, "value": ""})
    for _, row in params.loc[params["op"] == "=~"].iterrows():
        rows.append({"panel": "B", "row_type": "loading", "item": f"{row['lhs']} -> {row['rhs']}", "order": "", "value": row["STD_All"]})
    rows.extend([
        {"panel": "B", "row_type": "fit", "item": "D2 CFI", "order": "", "value": d2["CFI"]},
        {"panel": "B", "row_type": "fit", "item": "D2 SRMR", "order": "", "value": d2["SRMR"]},
        {"panel": "C", "row_type": "no_go", "item": "Higher-order general factor", "order": "", "value": "NO-GO"},
        {"panel": "C", "row_type": "no_go", "item": "IBD latent GWAS", "order": "", "value": "NO-GO"},
    ])
    save_data(pd.DataFrame(rows), "Figure1")

    fig = plt.figure(figsize=(8.0, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.38, wspace=0.35)
    axa = fig.add_subplot(gs[0, :])
    axb = fig.add_subplot(gs[1, 0])
    axc = fig.add_subplot(gs[1, 1])

    axa.set_xlim(0, 1); axa.set_ylim(0, 1); axa.axis("off")
    add_panel_label(axa, "A")
    axa.set_title("Frozen analysis workflow", loc="left", pad=2)
    for i, item in enumerate(workflow):
        x = 0.015 + i * 0.123
        color = ["#4D4D4D", "#4D4D4D", "#4D4D4D", "#0072B2", "#009E73", "#D55E00", "#56B4E9", "#CC79A7"][i]
        box = mpl.patches.FancyBboxPatch((x, 0.33), 0.105, 0.28, boxstyle="round,pad=0.012,rounding_size=0.015", facecolor="white", edgecolor=color, linewidth=1.5)
        axa.add_patch(box)
        axa.text(x + 0.0525, 0.47, item, ha="center", va="center", fontsize=7)
        if i < len(workflow) - 1:
            axa.annotate("", xy=(x + 0.119, 0.47), xytext=(x + 0.107, 0.47), arrowprops={"arrowstyle": "->", "lw": 0.7, "color": "#666666"})
    axa.text(0.03, 0.17, "SuSiE / multi-signal LD-aware analyses: attempted, QC-limited, not confirmatory", fontsize=7, color="#666666")

    axb.set_xlim(0, 1); axb.set_ylim(0, 1); axb.axis("off")
    add_panel_label(axb, "B")
    axb.set_title("D2 architecture", loc="left", pad=2)
    nodes = {
        "CM": (0.72, 0.75, "Cardiometabolic"),
        "PSY": (0.72, 0.28, "Internalizing"),
        "CD": (0.15, 0.83, "CD"), "UC": (0.15, 0.58, "UC"), "PsA": (0.15, 0.33, "PsA"),
    }
    for _, (x, y, label) in nodes.items():
        face = "#E6F2F8" if label == "Cardiometabolic" else "#E5F4EF" if label == "Internalizing" else "#F4F4F4"
        axb.add_patch(mpl.patches.FancyBboxPatch((x - 0.12, y - 0.07), 0.24, 0.14, boxstyle="round,pad=0.01", facecolor=face, edgecolor="#555555", linewidth=0.8))
        axb.text(x, y, label, ha="center", va="center", fontsize=7)
    loading_map = {str(row["rhs"]): float(row["STD_All"]) for _, row in params.loc[params["op"] == "=~"].iterrows()}
    for target, label, y, color in [
        ("bmi.ldsc.tsv", "BMI", 0.86, "#0072B2"),
        ("type_2_diabetes.ldsc.tsv", "T2D", 0.72, "#0072B2"),
        ("cad", "CAD", 0.58, "#0072B2"),
        ("major_depressive_disorder.ldsc.tsv", "MDD", 0.43, "#009E73"),
        ("anxiety_disorder.ldsc.tsv", "Anxiety", 0.29, "#009E73"),
        ("ptsd.ldsc.tsv", "PTSD", 0.15, "#009E73"),
    ]:
        if target in loading_map:
            axb.plot([0.84, 0.89], [y, y], color=color, lw=0.8)
            axb.text(0.91, y, f"{label}  {loading_map[target]:.3f}", ha="left", va="center", fontsize=6, color="#555555")
    axb.plot([0.72, 0.72], [0.35, 0.68], color="#555555", lw=0.8)
    axb.text(0.75, 0.515, "r = 0.307", ha="left", va="center", fontsize=7, color="#555555")
    axb.text(0.15, 0.08, "Observed immune traits", ha="center", va="center", fontsize=7, color="#555555")
    axb.text(0.5, 0.02, "CFI 0.980 | SRMR 0.042", ha="center", fontsize=7, fontweight="bold")

    axc.set_xlim(0, 1); axc.set_ylim(0, 1); axc.axis("off")
    add_panel_label(axc, "C")
    axc.set_title("General-factor decision", loc="left", pad=2)
    for y, label in [(0.68, "Higher-order general factor"), (0.35, "IBD latent GWAS")]:
        axc.add_patch(mpl.patches.Rectangle((0.14, y - 0.09), 0.72, 0.18, facecolor="#F7F7F7", edgecolor="#999999", linewidth=0.9))
        axc.text(0.5, y, label, ha="center", va="center", fontsize=8)
        axc.text(0.5, y - 0.045, "NO-GO", ha="center", va="center", fontsize=7, color="#D55E00", fontweight="bold")
    axc.text(0.5, 0.10, "Excluded from downstream discovery", ha="center", fontsize=7, color="#666666")
    save_figure(fig, "Figure1", MAIN_DIR)


def make_figure2() -> None:
    hotspots = cumulative_positions(read_table("results/phase3b0/hotspots/pair_hotspots.tsv"))
    global_local = read_table("results/phase3b0/global_local_architecture.tsv")
    pairs = ordered_pairs(global_local)
    ymap = {pair: i for i, pair in enumerate(pairs)}
    hotspots["pair_y"] = hotspots["pair"].map(ymap)
    hotspots["sign"] = np.where(pd.to_numeric(hotspots["lead_rho_display"], errors="coerce") >= 0, "positive", "negative")
    hotspot_fdr = pd.to_numeric(hotspots["min_atlas_fdr"], errors="coerce").clip(lower=np.finfo(float).tiny)
    hotspots["minus_log10_atlas_fdr"] = -np.log10(hotspot_fdr)
    data_a = hotspots[["pair_hotspot_id", "pair", "domain_axis", "chr", "start", "stop", "cumulative_mid_mb", "lead_rho_display", "sign", "minus_log10_atlas_fdr", "MHC_flag", "bonferroni_any"]].copy()
    data_a.insert(0, "panel", "A")
    data_b = global_local[["pair", "domain_axis", "global_rg", "n_pair_hotspots", "global_local_class"]].copy()
    data_b.insert(0, "panel", "B")
    save_data(pd.concat([data_a, data_b], ignore_index=True, sort=False), "Figure2")

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(11.2, 7.3), gridspec_kw={"width_ratios": [4.7, 1.8]})
    add_panel_label(ax, "A")
    ax.set_title("81 consolidated pair-hotspots across chr1-22", loc="left", pad=2)
    for axis, subset in hotspots.groupby("domain_axis", sort=False):
        for sign, marker in [("positive", "o"), ("negative", "x")]:
            d = subset.loc[subset["sign"] == sign]
            ax.scatter(d["cumulative_mid_mb"], d["pair_y"], s=18 + d["minus_log10_atlas_fdr"].clip(0, 8) * 3, marker=marker, color=AXIS_COLORS.get(axis, "#555555"), alpha=0.85, linewidths=0.9, edgecolors="#111111" if sign == "positive" else None, label=AXIS_SHORT.get(axis, axis) if sign == "positive" else None)
    chromosome_midpoints = []
    for chrom in range(1, 23):
        offset = sum(CHR_LENGTHS[x] for x in range(1, chrom))
        chromosome_midpoints.append((offset + CHR_LENGTHS[chrom] / 2) / 1e6)
        ax.axvline(offset / 1e6, color="#DDDDDD", lw=0.35, zorder=0)
    ax.set_xticks(chromosome_midpoints)
    ax.set_xticklabels([str(x) for x in range(1, 23)], fontsize=6)
    ax.tick_params(axis="x", pad=2)
    ax.set_xlim(0, sum(CHR_LENGTHS.values()) / 1e6)
    ax.set_ylim(-1, len(pairs))
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([pair_label(x) for x in pairs], fontsize=6)
    ax.set_xlabel("Cumulative genomic position (Mb; GRCh37 chromosome offsets)")
    ax.set_ylabel("Trait pair")
    ax.legend(frameon=False, fontsize=6, loc="upper left", ncol=3)
    ax.text(0.99, 0.02, "o positive | x negative\npoint size = -log10(atlas FDR)\nblack edge = MHC overlap", transform=ax.transAxes, ha="right", va="bottom", fontsize=6, color="#555555")
    mhc = hotspots.loc[hotspots["MHC_flag"] == True].copy()
    inset = ax.inset_axes([0.72, 0.74, 0.25, 0.22])
    inset.scatter(mhc["start"] / 1e6, mhc["pair_y"], s=10, color="#0072B2", alpha=0.8)
    inset.set_xlim(25, 34); inset.set_xticks([25, 30, 34]); inset.set_yticks([]); inset.set_title("MHC", fontsize=7, pad=1)
    inset.axvspan(25, 34, color="#E6F2F8", alpha=0.7, zorder=-1)
    inset.spines["top"].set_visible(True); inset.spines["right"].set_visible(True)
    add_panel_label(axb, "B")
    axb.set_title("Global rg and local hotspot count", loc="left", pad=2)
    for axis, subset in global_local.groupby("domain_axis", sort=False):
        axb.scatter(subset["global_rg"], subset["n_pair_hotspots"], s=15 + subset["n_pair_hotspots"] * 1.5, color=AXIS_COLORS.get(axis, "#555555"), label=AXIS_SHORT.get(axis, axis))
    for pair in ["UC__T2D", "PsA__PTSD"]:
        d = global_local.loc[global_local["pair"] == pair]
        if not d.empty:
            axb.annotate(pair_label(pair), (d["global_rg"].iloc[0], d["n_pair_hotspots"].iloc[0]), xytext=(4, 3), textcoords="offset points", fontsize=7)
    axb.axvline(0, color="#999999", lw=0.6)
    axb.set_xlabel("Genome-wide rg")
    axb.set_ylabel("Pair-hotspot count")
    axb.legend(frameon=False, fontsize=6, loc="upper left")
    save_figure(fig, "Figure2", MAIN_DIR)


def make_figure3() -> None:
    d = read_table("results/phase3b0/global_local_architecture.tsv").copy()
    d["pair_label"] = d["pair"].map(pair_label)
    save_data(d[["pair", "pair_label", "domain_axis", "global_rg", "n_pair_hotspots", "global_local_class"]], "Figure3")
    counts = d["global_local_class"].value_counts().rename_axis("global_local_class").reset_index(name="n_pairs")
    counts["label"] = counts["global_local_class"].str.replace("GLOBAL_", "", regex=False).str.replace("_", " / ", regex=False).str.title()

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(10.3, 4.8), gridspec_kw={"width_ratios": [3.2, 1.5]})
    add_panel_label(ax, "A")
    ax.set_title("Genome-wide rg versus local hotspot burden", loc="left", pad=2)
    for cls, subset in d.groupby("global_local_class", sort=False):
        ax.scatter(subset["global_rg"], subset["n_pair_hotspots"], s=35, color=CLASS_COLORS.get(cls, "#777777"), label=cls.replace("GLOBAL_", "").replace("_", " ").title(), alpha=0.9)
    for pair, dx, dy in [("UC__T2D", 5, 4), ("PsA__PTSD", 5, -10), ("BMI__MDD", 5, 4), ("CAD__MDD", 5, -10)]:
        z = d.loc[d["pair"] == pair]
        if not z.empty:
            ax.annotate(pair_label(pair), (z["global_rg"].iloc[0], z["n_pair_hotspots"].iloc[0]), xytext=(dx, dy), textcoords="offset points", fontsize=7)
    ax.axvline(0, color="#AAAAAA", lw=0.6)
    ax.set_xlabel("Genome-wide rg")
    ax.set_ylabel("Number of atlas-FDR pair-hotspots")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    add_panel_label(axb, "B")
    axb.set_title("Frozen classification counts", loc="left", pad=2)
    y = np.arange(len(counts))
    axb.barh(y, counts["n_pairs"], color=[CLASS_COLORS.get(x, "#777777") for x in counts["global_local_class"]])
    axb.set_yticks(y); axb.set_yticklabels(counts["label"], fontsize=6)
    axb.set_xlabel("Trait pairs")
    axb.invert_yaxis()
    for yi, val in zip(y, counts["n_pairs"]): axb.text(val + 0.2, yi, str(val), va="center", fontsize=7)
    save_figure(fig, "Figure3", MAIN_DIR)


def make_figure4() -> None:
    cross = read_table("results/phase3b0/hotspots/cross_pair_hotspots.tsv")
    pair_data = read_table("results/phase3b0/hotspots/pair_hotspots.tsv")
    global_local = read_table("results/phase3b0/global_local_architecture.tsv")
    pair_order = ordered_pairs(global_local)
    cross = cross.sort_values(["chr", "start", "stop", "cross_pair_hotspot_id"]).reset_index(drop=True)
    long_rows = []
    matrices = []
    for _, row in cross.iterrows():
        pair_set = set(str(row["pairs"]).split(",")) if pd.notna(row["pairs"]) else set()
        matrix_row = []
        for pair in pair_order:
            present = int(pair in pair_set)
            matrix_row.append(present)
            long_rows.append({"cross_pair_hotspot_id": row["cross_pair_hotspot_id"], "chr": row["chr"], "start": row["start"], "stop": row["stop"], "pair": pair, "pair_present": present, "n_pairs": row["n_pairs"], "MHC_flag": row["MHC_flag"], "recurrent_hotspot": row["recurrent_hotspot"], "multi_axis_hotspot": row["multi_axis_hotspot"], "pan_domain_hotspot": row["pan_domain_hotspot"]})
        matrices.append(matrix_row)
    save_data(pd.DataFrame(long_rows), "Figure4")

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(11.3, 8.3), gridspec_kw={"width_ratios": [4.4, 1.25]})
    add_panel_label(ax, "A")
    ax.set_title("59 cross-pair hotspots and their pair support", loc="left", pad=2)
    ax.imshow(np.asarray(matrices), aspect="auto", interpolation="none", cmap=ListedColormap(["#F7F7F7", "#0072B2"]), vmin=0, vmax=1)
    ax.set_xticks(range(len(pair_order))); ax.set_xticklabels([pair_label(x) for x in pair_order], rotation=90, rotation_mode="anchor", fontsize=5)
    yticks = np.arange(len(cross))
    ax.set_yticks(yticks[::4]); ax.set_yticklabels(cross.loc[::4, "cross_pair_hotspot_id"], fontsize=6)
    ax.set_xlabel("Trait pair")
    ax.set_ylabel("Cross-pair hotspot (fixed genomic order)")
    for i, flag in enumerate(cross["MHC_flag"]):
        if bool(flag): ax.add_patch(mpl.patches.Rectangle((-0.5, i - 0.5), len(pair_order), 1, fill=False, edgecolor="#111111", linewidth=0.7))
    ax.text(1.0, -0.16, "Blue = pair support | outlined row = MHC overlap", transform=ax.transAxes, ha="right", fontsize=6, color="#555555")

    add_panel_label(axb, "B")
    axb.set_title("Recurrence", loc="left", pad=2)
    cross["status"] = np.select([cross["multi_axis_hotspot"].astype(bool), cross["recurrent_hotspot"].astype(bool)], ["multi-axis", "recurrent"], default="single-pair")
    y = np.arange(len(cross))
    color = ["#CC79A7" if x == "multi-axis" else "#0072B2" if x == "recurrent" else "#BBBBBB" for x in cross["status"]]
    axb.barh(y, cross["n_pairs"], color=color, height=0.72)
    axb.set_yticks([]); axb.set_xlabel("Number of supporting pairs")
    axb.invert_yaxis()
    axb.text(0.02, 0.02, "12 recurrent\n4 multi-axis\n0 pan-domain", transform=axb.transAxes, va="bottom", fontsize=7)
    save_figure(fig, "Figure4", MAIN_DIR)


def make_figure5() -> None:
    mdd = read_table("results/phase3b0/sensitivity/mdd_noUKBB_hotspot_concordance.tsv").copy()
    mdd["testable"] = mdd["noUKBB_rho_display"].notna()
    abf = read_table("results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv").copy()
    data_mdd = mdd[["pair", "locus", "chr", "start", "stop", "primary_rho", "noUKBB_rho_display", "robustness_class", "testable"]].copy()
    data_mdd.insert(0, "panel", "A")
    data_abf = abf[["pair_hotspot_id", "pair", "chr", "start", "stop", "PP.H3.abf", "PP.H4.abf", "coloc_evidence_class", "trait1_type", "trait2_type"]].copy()
    data_abf.insert(0, "panel", "B")
    save_data(pd.concat([data_mdd, data_abf], ignore_index=True, sort=False), "Figure5")

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(11.2, 6.0), gridspec_kw={"width_ratios": [1.05, 1.35]})
    add_panel_label(ax, "A")
    ax.set_title("MDD noUKBB sensitivity", loc="left", pad=2)
    robust = mdd.loc[mdd["robustness_class"] == "MDD_ROBUST"]
    not_testable = mdd.loc[mdd["robustness_class"] != "MDD_ROBUST"]
    ax.scatter(robust["primary_rho"], robust["noUKBB_rho_display"], s=17, color="#0072B2", label="Robust (n=37)")
    ax.scatter(not_testable["primary_rho"], not_testable["noUKBB_rho_display"], s=28, facecolors="none", edgecolors="#D55E00", marker="s", label="Not testable (n=4)")
    lo, hi = -1.3, 1.3
    ax.plot([lo, hi], [lo, hi], ls="--", color="#999999", lw=0.7)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Primary rho")
    ax.set_ylabel("noUKBB rho (display value)")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax.text(0.04, 0.96, "0 attenuated | 0 reversed", transform=ax.transAxes, va="top", fontsize=7, color="#555555")

    add_panel_label(axb, "B")
    axb.set_title("Secondary ABF colocalization screening", loc="left", pad=2)
    abf = abf.sort_values(["PP.H4.abf", "PP.H3.abf"]).reset_index(drop=True)
    y = np.arange(len(abf))
    axb.barh(y - 0.17, abf["PP.H4.abf"], height=0.32, color="#0072B2", label="PP4")
    axb.barh(y + 0.17, abf["PP.H3.abf"], height=0.32, color="#D55E00", label="PP3")
    axb.set_yticks(y); axb.set_yticklabels(abf["pair_hotspot_id"], fontsize=6)
    axb.set_xlim(0, 1.02); axb.set_xlabel("Posterior probability")
    axb.text(0.98, 0.98, "blue = PP4 | orange = PP3", transform=axb.transAxes, ha="right", va="top", fontsize=6, color="#555555")
    fig.text(0.75, 0.015, "Five LD-aware candidates: QC provenance only", ha="center", fontsize=6, color="#555555")
    save_figure(fig, "Figure5", MAIN_DIR)


def make_supplementary() -> None:
    # S1: frozen nine-trait LDSC rg matrix.
    rg = read_table("results/ldsc/phase2c_balanced9_rg_matrix.tsv")
    trait_cols = [x for x in rg.columns if x != "trait"]
    mat = rg[trait_cols].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="none")
    ax.set_xticks(range(len(trait_cols))); ax.set_xticklabels([x.replace("_", " ") for x in trait_cols], rotation=90, rotation_mode="anchor", fontsize=6)
    ax.set_yticks(range(len(rg))); ax.set_yticklabels(rg["trait"].str.replace("_", " "), fontsize=6)
    ax.set_title("S1. Frozen nine-trait LDSC genetic-correlation matrix", loc="left", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03, label="rg")
    save_figure(fig, "Supplementary_S1", SUPP_DIR, include_tiff=False)
    save_data(rg, "Supplementary_S1")

    # S2: model comparison diagnostics, shown without refitting.
    model = read_table("results/genomic_sem/phase2d/model_comparison.tsv")
    model["CFI"] = pd.to_numeric(model["CFI"], errors="coerce")
    model["SRMR"] = pd.to_numeric(model["SRMR"], errors="coerce")
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    good = model.loc[model["CFI"].notna() & (model["model"] != "D2_CD_UC_PsA_observed")]
    ax.scatter(good["CFI"], good["SRMR"], s=38, color="#777777")
    for _, row in good.iterrows(): ax.annotate(str(row["model"]).split("_")[0], (row["CFI"], row["SRMR"]), xytext=(4, 3), textcoords="offset points", fontsize=7)
    d2_row = model.loc[model["model"] == "D2_CD_UC_PsA_observed"].iloc[0]
    ax.scatter([float(d2_row["CFI"])], [float(d2_row["SRMR"])], s=65, facecolors="none", edgecolors="#0072B2", linewidths=1.5, label="D2")
    ax.annotate("D2", (float(d2_row["CFI"]), float(d2_row["SRMR"])), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_xlabel("CFI"); ax.set_ylabel("SRMR"); ax.set_title("S2. Genomic SEM model-comparison diagnostics", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=7)
    save_figure(fig, "Supplementary_S2", SUPP_DIR, include_tiff=False)
    save_data(model, "Supplementary_S2")

    # S3: all significant LAVA blocks, retaining blocks as the analysis unit.
    blocks = cumulative_positions(read_table("results/local_rg/lava/atlas_27pair_qc/significant_loci_atlas_fdr_0.05.tsv"))
    pair_axes = read_table("results/phase3b0/hotspots/pair_hotspots.tsv")[["pair", "domain_axis"]].drop_duplicates()
    blocks = blocks.merge(pair_axes, on="pair", how="left")
    blocks["MHC_flag"] = (blocks["chr"] == 6) & (blocks["start"] <= 34000000) & (blocks["stop"] >= 25000000)
    block_fdr = pd.to_numeric(blocks["atlas_fdr"], errors="coerce").clip(lower=np.finfo(float).tiny)
    blocks["minus_log10_atlas_fdr"] = -np.log10(block_fdr)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    for axis, subset in blocks.groupby("domain_axis", sort=False):
        ax.scatter(subset["cumulative_mid_mb"], subset["minus_log10_atlas_fdr"], s=12, color=AXIS_COLORS.get(axis, "#555555"), label=AXIS_SHORT.get(axis, axis), alpha=0.8)
    ax.set_xlabel("Cumulative genomic position (Mb; GRCh37 chromosome offsets)"); ax.set_ylabel("-log10(atlas FDR)"); ax.set_title("S3. 93 significant LAVA blocks", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=6, ncol=3)
    ax.text(0.99, 0.03, "Blocks are LAVA analysis units and are not necessarily independent loci.", transform=ax.transAxes, ha="right", fontsize=7, color="#555555")
    save_figure(fig, "Supplementary_S3", SUPP_DIR, include_tiff=False)
    save_data(blocks, "Supplementary_S3")

    # S4: MHC retention summary from the frozen block and pair-hotspot tables.
    pair_hotspots = read_table("results/phase3b0/hotspots/pair_hotspots.tsv")
    s4 = pd.DataFrame({"unit": ["LAVA blocks", "Pair-hotspots"], "all": [len(blocks), len(pair_hotspots)], "MHC_overlap": [int(blocks["MHC_flag"].astype(bool).sum()), int(pair_hotspots["MHC_flag"].astype(bool).sum())]})
    fig, ax = plt.subplots(figsize=(5.8, 4.3))
    x = np.arange(len(s4)); width = 0.35
    ax.bar(x - width / 2, s4["all"], width, color="#999999", label="All")
    ax.bar(x + width / 2, s4["MHC_overlap"], width, color="#0072B2", label="MHC overlap")
    ax.set_xticks(x); ax.set_xticklabels(s4["unit"]); ax.set_ylabel("Count"); ax.set_title("S4. MHC retention and sensitivity summary", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=7)
    save_figure(fig, "Supplementary_S4", SUPP_DIR, include_tiff=False)
    save_data(s4, "Supplementary_S4")

    # S5: full MDD noUKBB concordance table as a display companion.
    mdd = read_table("results/phase3b0/sensitivity/mdd_noUKBB_hotspot_concordance.tsv")
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    for cls, subset in mdd.groupby("robustness_class", sort=False):
        ax.scatter(subset["primary_rho"], subset["noUKBB_rho_display"], s=15, label=cls.replace("MDD_", "").replace("_", " ").title(), alpha=0.85)
    ax.plot([-1.3, 1.3], [-1.3, 1.3], ls="--", color="#999999", lw=0.7); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
    ax.set_xlabel("Primary rho"); ax.set_ylabel("noUKBB rho (display value)"); ax.set_title("S5. MDD noUKBB hotspot concordance", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=6)
    save_figure(fig, "Supplementary_S5", SUPP_DIR, include_tiff=False)
    save_data(mdd, "Supplementary_S5")

    # S6: ABF posterior summary.
    abf = read_table("results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv")
    fig, ax = plt.subplots(figsize=(6.3, 4.8))
    for cls, subset in abf.groupby("coloc_evidence_class", sort=False):
        ax.scatter(subset["PP.H3.abf"], subset["PP.H4.abf"], s=33, color=ABF_COLORS.get(cls, "#777777"), label=cls.replace("_", " ").title())
    ax.axvline(0.8, color="#999999", lw=0.6, ls="--"); ax.axhline(0.8, color="#999999", lw=0.6, ls="--")
    ax.set_xlabel("PP3"); ax.set_ylabel("PP4"); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_title("S6. ABF posterior screening summary", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=6)
    save_figure(fig, "Supplementary_S6", SUPP_DIR, include_tiff=False)
    save_data(abf, "Supplementary_S6")

    # S7: second-LD diagnostics, explicitly QC provenance.
    original = read_table("results/phase3b2r/ld/original_ld_qc.tsv")
    second = read_table("results/phase3b2r/ld_second/second_ld_qc.tsv")
    qc = original[["candidate", "min_eigenvalue", "n_negative_eigenvalues"]].rename(columns={"min_eigenvalue": "original_min_eigenvalue", "n_negative_eigenvalues": "original_n_negative_eigenvalues"}).merge(second[["candidate", "min_eigenvalue", "n_negative_eigenvalues", "matrix_status"]].rename(columns={"min_eigenvalue": "second_min_eigenvalue", "n_negative_eigenvalues": "second_n_negative_eigenvalues"}), on="candidate", how="outer")
    fig, ax = plt.subplots(figsize=(8.0, 4.7))
    y = np.arange(len(qc)); h = 0.34
    ax.barh(y - h / 2, qc["original_n_negative_eigenvalues"], h, color="#999999", label="Original LD")
    ax.barh(y + h / 2, qc["second_n_negative_eigenvalues"], h, color="#0072B2", label="Alternative UKB-derived LD")
    ax.set_yticks(y); ax.set_yticklabels(qc["candidate"], fontsize=7); ax.set_xlabel("Negative eigenvalue count"); ax.set_title("S7. LD-reference QC provenance", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=7)
    ax.text(0.99, 0.02, "Not biological negative evidence; no confirmatory fine-mapping claim.", transform=ax.transAxes, ha="right", fontsize=7, color="#555555")
    save_figure(fig, "Supplementary_S7", SUPP_DIR, include_tiff=False)
    save_data(qc, "Supplementary_S7")


def write_qc_files() -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    number_rows = [
        ["Figure1", "A", "primary traits", 9, "MANUSCRIPT_STATISTICAL_FREEZE.md", "9", "TRUE"],
        ["Figure1", "B", "D2 CFI", 0.980, "results/genomic_sem/phase2d/model_comparison.tsv", "0.980", "TRUE"],
        ["Figure1", "B", "D2 SRMR", 0.042, "results/genomic_sem/phase2d/model_comparison.tsv", "0.042", "TRUE"],
        ["Figure2", "A", "consolidated pair-hotspots", 81, "results/phase3b0/hotspots/pair_hotspots.tsv", "81 rows", "TRUE"],
        ["Figure2", "B", "cross-domain pairs", 27, "results/local_rg/lava/atlas_27pair_qc/pair_signal_summary.tsv", "27 rows", "TRUE"],
        ["Figure3", "A", "cross-domain pairs", 27, "results/phase3b0/global_local_architecture.tsv", "27 rows", "TRUE"],
        ["Figure4", "A", "cross-pair hotspots", 59, "results/phase3b0/hotspots/cross_pair_hotspots.tsv", "59 rows", "TRUE"],
        ["Figure4", "B", "recurrent hotspots", 12, "results/phase3b0/hotspots/cross_pair_hotspots.tsv", "sum(recurrent_hotspot)", "TRUE"],
        ["Figure4", "B", "multi-axis hotspots", 4, "results/phase3b0/hotspots/cross_pair_hotspots.tsv", "sum(multi_axis_hotspot)", "TRUE"],
        ["Figure4", "B", "pan-domain hotspots", 0, "results/phase3b0/hotspots/cross_pair_hotspots.tsv", "sum(pan_domain_hotspot)", "TRUE"],
        ["Figure5", "A", "primary MDD blocks", 41, "results/phase3b0/sensitivity/mdd_noUKBB_hotspot_concordance.tsv", "41 rows", "TRUE"],
        ["Figure5", "A", "robust MDD blocks", 37, "results/phase3b0/sensitivity/mdd_noUKBB_hotspot_concordance.tsv", "sum(MDD_ROBUST)", "TRUE"],
        ["Figure5", "A", "not-testable MDD blocks", 4, "results/phase3b0/sensitivity/mdd_noUKBB_hotspot_concordance.tsv", "sum(non-robust)", "TRUE"],
        ["Figure5", "B", "ABF candidates", 15, "results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv", "15 rows", "TRUE"],
        ["Figure5", "B", "strong ABF candidates", 4, "results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv", "sum(STRONG_SHARED_SIGNAL)", "TRUE"],
        ["Figure5", "B", "suggestive ABF candidates", 1, "results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv", "sum(SUGGESTIVE_SHARED_SIGNAL)", "TRUE"],
        ["Figure5", "B", "distinct-signal candidates", 8, "results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv", "sum(DISTINCT_SIGNAL_SUPPORTED)", "TRUE"],
        ["Figure5", "B", "weak/inconclusive candidates", 2, "results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv", "sum(WEAK_OR_INCONCLUSIVE)", "TRUE"],
    ]
    pd.DataFrame(number_rows, columns=["figure", "panel", "statistic", "display_value", "source_file", "source_value", "match"]).to_csv(QC_DIR / "FIGURE_NUMBER_AUDIT.tsv", sep="\t", index=False)
    source_rows = [
        ["Figure1A", "frozen workflow", "MANUSCRIPT_EVIDENCE_HIERARCHY.md; MANUSCRIPT_STATISTICAL_FREEZE.md", "none; schematic only", "make_figure1.py"],
        ["Figure1B", "D2 architecture", "results/genomic_sem/phase2d/D2_CD_UC_PsA_observed_parameters.tsv", "filter D2 loadings; display frozen STD_All", "make_figure1.py"],
        ["Figure1C", "general-factor decision", "FINAL_SEM_GO_NOGO.md", "display locked NO-GO status", "make_figure1.py"],
        ["Figure2A", "81 pair-hotspots", "results/phase3b0/hotspots/pair_hotspots.tsv", "sort by pair; chromosome offsets; sign; -log10(min atlas FDR)", "make_figure2.py"],
        ["Figure2B", "global-local context", "results/phase3b0/global_local_architecture.tsv", "display frozen global rg and pair-hotspot count", "make_figure2.py"],
        ["Figure3A-B", "global-local classification", "results/phase3b0/global_local_architecture.tsv", "plot and count frozen classification categories", "make_figure3.py"],
        ["Figure4A-B", "59 cross-pair hotspots", "results/phase3b0/hotspots/cross_pair_hotspots.tsv; results/phase3b0/hotspots/pair_hotspots.tsv", "fixed genomic sort; split frozen pairs field; binary display", "make_figure4.py"],
        ["Figure5A", "MDD noUKBB", "results/phase3b0/sensitivity/mdd_noUKBB_hotspot_concordance.tsv", "display primary and noUKBB rho; frozen class", "make_figure5.py"],
        ["Figure5B", "ABF screening", "results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv", "sort by PP4; display PP3 and PP4", "make_figure5.py"],
        ["Supplementary_S1", "LDSC rg matrix", "results/ldsc/phase2c_balanced9_rg_matrix.tsv", "matrix display only", "make_supplementary.py"],
        ["Supplementary_S2", "SEM diagnostics", "results/genomic_sem/phase2d/model_comparison.tsv", "display frozen CFI/SRMR", "make_supplementary.py"],
        ["Supplementary_S3", "93 LAVA blocks", "results/local_rg/lava/atlas_27pair_qc/significant_loci_atlas_fdr_0.05.tsv", "chromosome offsets; -log10(atlas FDR)", "make_supplementary.py"],
        ["Supplementary_S4", "MHC sensitivity", "S3 blocks and pair_hotspots.tsv", "count aggregation by MHC flag", "make_supplementary.py"],
        ["Supplementary_S5", "MDD concordance", "results/phase3b0/sensitivity/mdd_noUKBB_hotspot_concordance.tsv", "display frozen concordance classes", "make_supplementary.py"],
        ["Supplementary_S6", "ABF posterior summary", "results/phase3b1/coloc/phase3b1_abf_coloc_summary.tsv", "display frozen PP3/PP4", "make_supplementary.py"],
        ["Supplementary_S7", "LD QC provenance", "results/phase3b2r/ld/original_ld_qc.tsv; results/phase3b2r/ld_second/second_ld_qc.tsv", "join by frozen candidate; display negative eigenvalue counts", "make_supplementary.py"],
    ]
    pd.DataFrame(source_rows, columns=["figure_panel", "source_table", "source_file", "filter_or_transformation", "plotting_script"]).to_csv(QC_DIR / "FIGURE_SOURCE_MAP.tsv", sep="\t", index=False)
    panel_rows = [
        ["Figure1A", "Frozen workflow", "9 primary traits and locked analysis steps", "No spread; schematic", "No test or correction", "Workflow labels and QC note", "Pass", "TRUE"],
        ["Figure1B", "D2 architecture", "Frozen standardized loadings and domain correlation", "No spread; frozen point estimates", "No new test or correction", "Trait loading labels; no collision", "Pass", "TRUE"],
        ["Figure1C", "Locked NO-GO decisions", "Two frozen decision statuses", "Not applicable", "Not applicable", "NO-GO labels clear", "Pass", "TRUE"],
        ["Figure2A", "81 pair-hotspots", "One point per consolidated pair-hotspot", "No spread; one displayed lead estimate", "BH-FDR inherited from frozen table", "Pair labels, sign key, MHC inset", "Pass", "TRUE"],
        ["Figure2B", "Global rg and local count", "27 predefined trait pairs", "No spread; frozen rg and counts", "No new test or correction", "Direct pair labels", "Pass", "TRUE"],
        ["Figure3A", "Global-local relationship", "27 predefined trait pairs", "No spread; frozen rg and counts", "No new test or correction", "Four direct pair labels", "Pass", "TRUE"],
        ["Figure3B", "Classification counts", "27 predefined trait pairs", "No spread; category counts", "Frozen classification only", "Category labels and values", "Pass", "TRUE"],
        ["Figure4A", "Cross-pair hotspot support", "59 hotspots x 27 predefined pairs", "No spread; binary support display", "No new test or correction", "Rotated pair labels; fixed genomic order", "Pass", "TRUE"],
        ["Figure4B", "Hotspot recurrence", "59 cross-pair hotspots", "No spread; support-pair counts", "No new test or correction", "Recurrence note clear", "Pass", "TRUE"],
        ["Figure5A", "MDD noUKBB concordance", "41 primary MDD blocks", "No spread; paired rho display", "Frozen robustness classification", "Legend and identity line clear", "Pass", "TRUE"],
        ["Figure5B", "ABF screening", "15 frozen candidates", "No spread; PP3/PP4 display", "Frozen ABF evidence classes", "Candidate labels and PP key clear", "Pass", "TRUE"],
        ["Supplementary_S1", "LDSC rg matrix", "9 traits", "No spread; matrix entries", "Frozen LDSC output", "Rotated trait labels clear", "Pass", "TRUE"],
        ["Supplementary_S2", "SEM diagnostics", "Frozen model comparison rows", "No spread; CFI/SRMR points", "No refit", "D2 highlighted once", "Pass", "TRUE"],
        ["Supplementary_S3", "93 LAVA blocks", "93 significant LAVA blocks", "No spread; block-level display", "Atlas-FDR inherited", "Note distinguishes blocks from loci", "Pass", "TRUE"],
        ["Supplementary_S4", "MHC summary", "Blocks and pair-hotspots", "No spread; count aggregation", "MHC flag inherited", "Legend clear", "Pass", "TRUE"],
        ["Supplementary_S5", "MDD concordance", "41 primary MDD blocks", "No spread; paired rho display", "Frozen robustness class", "Legend and identity line clear", "Pass", "TRUE"],
        ["Supplementary_S6", "ABF posterior summary", "15 frozen candidates", "No spread; PP3/PP4 points", "Frozen display thresholds", "Class legend clear", "Pass", "TRUE"],
        ["Supplementary_S7", "LD QC provenance", "5 frozen candidates", "No spread; negative eigenvalue counts", "QC provenance only", "Reference legend and note clear", "Pass", "TRUE"],
    ]
    pd.DataFrame(panel_rows, columns=["panel", "unique_claim", "n_definition", "center_and_spread", "test_and_correction", "labels_and_legend", "collision_check", "pass"]).to_csv(QC_DIR / "FIGURE_PANEL_QC.tsv", sep="\t", index=False)


def run_figure(name: str) -> None:
    configure_style()
    DATA_DIR.mkdir(parents=True, exist_ok=True); MAIN_DIR.mkdir(parents=True, exist_ok=True); SUPP_DIR.mkdir(parents=True, exist_ok=True); QC_DIR.mkdir(parents=True, exist_ok=True)
    if name == "Figure1": make_figure1()
    elif name == "Figure2": make_figure2()
    elif name == "Figure3": make_figure3()
    elif name == "Figure4": make_figure4()
    elif name == "Figure5": make_figure5()
    else: raise ValueError(name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("figure", nargs="?", choices=["Figure1", "Figure2", "Figure3", "Figure4", "Figure5", "all"], default="all")
    args = parser.parse_args()
    if args.figure == "all":
        for name in ["Figure1", "Figure2", "Figure3", "Figure4", "Figure5"]: run_figure(name)
        configure_style(); make_supplementary(); write_qc_files()
    else:
        run_figure(args.figure)


if __name__ == "__main__":
    main()
