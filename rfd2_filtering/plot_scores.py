"""
Renders a histogram per numeric column in scores.csv into a timestamped
subfolder under rfd2_filtering/plots/, so you can pull the whole folder
down and skim distributions without a Python session.

Usage:
    python -m rfd2_filtering.plot_scores --scores-path rfd2_filtering/scores.csv
"""


import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
# PROJECT_DIR = SCRIPT_DIR.parent


def plot_column(df: pd.DataFrame, column: str, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    data = df[column].dropna()
    ax.hist(data, bins="auto")  # numpy/matplotlib's standard auto bin-width rule (max of Sturges and Freedman-Diaconis)
    ax.set_xlabel(column)
    ax.set_ylabel("count")
    ax.set_title(f"Distribution of {column} (n={len(data)})")
    return ax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores-path", type=Path, default=SCRIPT_DIR / "scores.csv")
    ap.add_argument("--plots-dir", type=Path, default=SCRIPT_DIR / "plots")
    args = ap.parse_args()

    df = pd.read_csv(args.scores_path)

    out_dir = args.plots_dir / datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir.mkdir(parents=True, exist_ok=True)

    skip_cols = {"design_path", "config_id"}
    numeric_cols = [c for c in df.columns if c not in skip_cols and pd.api.types.is_numeric_dtype(df[c])]

    n_written = 0
    for col in numeric_cols:
        if df[col].dropna().empty:
            print(f"Skipping {col}: no data")
            continue
        ax = plot_column(df, col)
        ax.figure.savefig(out_dir / f"{col}.png")
        ax.figure.clf()
        n_written += 1

    print(f"Wrote {n_written} plots to {out_dir}")


if __name__ == "__main__":
    main()

