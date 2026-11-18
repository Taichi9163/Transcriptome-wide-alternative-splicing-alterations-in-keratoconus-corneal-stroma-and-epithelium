"""Select Reactome top-ten pathways and optionally export the study bubble plots."""

import argparse
import re
import textwrap

from utils_paths import add_work_dir_argument

GENE_COLUMN = "Spliced genes (number of interactions affecting the pathway)"


def genes(value):
    """Return distinct gene symbols after removing interaction-count suffixes."""
    import pandas as pd

    if pd.isna(value):
        return []
    return sorted({re.sub(r"\s*\([^)]*\)", "", v).strip() for v in value.split(",") if v.strip()})


def select_top_pathways(frame):
    """Apply the inclusive BH cutoff, then rank by score and stable pathway ID."""
    import numpy as np

    selected = (
        frame.loc[frame["adj p_value"] <= 0.05]
        .sort_values(["Nease score", "Pathway ID"], ascending=[False, True])
        .head(10)
        .copy()
    )
    assert len(selected) == 10 and selected["Pathway ID"].is_unique
    assert selected["Nease score"].is_monotonic_decreasing
    selected["Splicing-affected genes"] = selected[GENE_COLUMN].map(lambda value: len(genes(value)))
    assert (selected["Splicing-affected genes"] > 0).all() and (selected["adj p_value"] > 0).all()
    selected["-log10 adjusted P"] = -np.log10(selected["adj p_value"])
    selected.insert(0, "Rank", range(1, 11))
    return selected


def plot_pathways(p, run, epithelial, large_text):
    """Render the original layout; called only when figure export is requested."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.ticker import MultipleLocator, MaxNLocator

    plt.rcParams.update(
        {"font.family": "Arial", "font.size": 12, "pdf.fonttype": 42, "ps.fonttype": 42}
    )
    (run / "figures").mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8.5, 6.5 if large_text else 6.0), facecolor="white")
    ax = fig.add_axes([0.515, 0.105, 0.25, 0.84])
    cmap = LinearSegmentedColormap.from_list(
        "nease", ["#FFFFFF", "#CC79A7" if epithelial else "#FAA62B"]
    )
    norm = Normalize(0, float(p["-log10 adjusted P"].max()))
    scale = 40 if epithelial else 16
    dots = ax.scatter(
        p["Nease score"],
        range(10),
        s=p["Splicing-affected genes"] * scale,
        c=p["-log10 adjusted P"],
        cmap=cmap,
        norm=norm,
        edgecolors="white",
        linewidths=0.7,
        zorder=3,
    )
    labels = [
        "\n".join(
            textwrap.wrap(
                v, 34 if large_text else 43, break_long_words=False, break_on_hyphens=False
            )
        )
        for v in p["Pathway name"]
    ]
    ax.set_yticks(range(10), labels)
    for label in ax.get_yticklabels():
        label.set_linespacing(1.03)
    ax.set_ylim(9.6, -0.6)
    lo, hi = p["Nease score"].min(), p["Nease score"].max()
    span = hi - lo
    ax.set_xlim(lo - span * 0.15, hi + span * 0.20)
    ax.set_xlabel("NEASE score", fontsize=18 if large_text else 14, labelpad=8)
    ax.tick_params(axis="both", labelsize=16 if large_text else 12, length=3, pad=5)
    ax.xaxis.set_major_locator(MultipleLocator(2 if epithelial else 10))
    ax.xaxis.set_minor_locator(MultipleLocator(1 if epithelial else 5))
    ax.grid(True, axis="both", which="major", linestyle=(0, (1, 3)), linewidth=0.7, color="#999999")
    ax.grid(True, axis="x", which="minor", linestyle=(0, (1, 3)), linewidth=0.5, color="#BBBBBB")
    for spine in ax.spines.values():
        spine.set_color("#555555")
        spine.set_linewidth(0.8)
    sizes = [2, 7, 14] if epithelial else [10, 20, 40]
    handles = [
        ax.scatter([], [], s=value * scale, color="black", edgecolors="none") for value in sizes
    ]
    fig.legend(
        handles,
        [str(v) for v in sizes],
        title="Splicing-affected\ngenes",
        loc="upper left",
        bbox_to_anchor=(0.79, 0.94),
        frameon=False,
        fontsize=14 if large_text else 11,
        title_fontsize=15 if large_text else 12,
        labelspacing=1.05,
        handletextpad=0.9,
        borderpad=0,
    )
    colorax = fig.add_axes([0.80, 0.19, 0.025, 0.29])
    cb = fig.colorbar(dots, cax=colorax)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=14 if large_text else 11, length=0, pad=6)
    cb.locator = MaxNLocator(nbins=5)
    cb.update_ticks()
    fig.text(
        0.80, 0.51, r"$-\mathrm{Log}_{10}$(P-adjust)", fontsize=15 if large_text else 12, ha="left"
    )
    for ext in ["png", "pdf"]:
        suffix = "_large_text" if large_text else ""
        target = run / "figures" / f"NEASE_Reactome_top10_updated{suffix}.{ext}"
        fig.savefig(target, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.08)
        print(target)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_work_dir_argument(parser)
    parser.add_argument(
        "--tables-only",
        action="store_true",
        help="Export source tables without opening graphics devices",
    )
    parser.add_argument("--large-text", action="store_true")
    args = parser.parse_args()
    import pandas as pd

    for cohort in ["stroma", "epithelium"]:
        run = args.work_dir / "results" / cohort / "nease"
        frame = pd.read_csv(run / "results/Reactome_Python.tsv", sep="\t")
        selected = select_top_pathways(frame)
        selected.to_csv(run / "results/Reactome_top10_for_figure.tsv", sep="\t", index=False)
        if not args.tables_only:
            plot_pathways(selected, run, cohort == "epithelium", args.large_text)


if __name__ == "__main__":
    main()
