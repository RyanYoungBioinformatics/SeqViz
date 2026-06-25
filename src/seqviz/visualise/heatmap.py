from pathlib import Path

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_matrix(
    matrix: np.ndarray,
    seq1: str,
    seq2: str,
    title: str = "Alignment Scoring Matrix",
    output_path: str | None = None,
) -> Path | None:
    """
    Render the scoring matrix as a heatmap.
    Saves to file if output_path is given, otherwise opens an interactive window.
    Returns the resolved absolute path if saved, otherwise None.
    """
    rows, cols = matrix.shape

    fig_width  = max(8, cols * 0.8)
    fig_height = max(6, rows * 0.8)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        matrix,
        xticklabels=["-"] + list(seq2),
        yticklabels=["-"] + list(seq1),
        annot=True,
        fmt="d",
        cmap="coolwarm",
        linewidths=0.5,
        linecolor="white",
        ax=ax,
    )

    ax.set_title(title, fontsize=14, pad=12)
    ax.set_xlabel("Sequence 2", fontsize=11)
    ax.set_ylabel("Sequence 1", fontsize=11)
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)

    plt.tight_layout()

    if output_path:
        resolved = Path(output_path).resolve()
        resolved.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(resolved, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return resolved       # ← this is what cli.py prints as the save location

    else:
        plt.show()
        return None
