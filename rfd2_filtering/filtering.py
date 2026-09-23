"""
Apply threshold filters to rfd2_filtering/scores.csv and write the
surviving designs to passing.csv.

Look at the histograms in rfd2_filtering/plots/<timestamp>/ (from
plot_scores.py) first to pick sensible bounds.

Usage:
    # list which columns you can filter on, and their observed min/max
    python -m rfd2_filtering.filtering --list-bounds

    # use the bounds hardcoded in BOUNDS below, write passing.csv
    python -m rfd2_filtering.filtering

    # try bounds ad-hoc without touching BOUNDS or writing passing.csv
    python -m rfd2_filtering.filtering --dry-run --bound n_clashes 0 0 --bound rg_A 10 25

    # same, but also write passing.csv with these bounds instead of BOUNDS
    python -m rfd2_filtering.filtering --bound n_clashes 0 0 --bound rg_A 10 25
"""

import argparse
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

# (min, max), inclusive. Used only if --bound isn't passed on the command line.
BOUNDS = {
    "n_clashes": (0, 0),
    "n_backbone_breaks": (0, 0),
}


def load_scores(scores_path: Path = SCRIPT_DIR / "scores.csv") -> pd.DataFrame:
    return pd.read_csv(scores_path, dtype={"config_id": str, "design_path": str})


def apply_filters(df: pd.DataFrame, bounds: dict) -> pd.DataFrame:
    """
    NaNs in a filtered column are dropped (treated as a fail, since NaN
    usually means the metric couldn't be computed for that design -- e.g.
    com_shift_A/B missing for a chain with no motif residues, or frac_loop
    missing because DSSP was skipped after a clash/backbone fail).
    """
    mask = pd.Series(True, index=df.index)
    for column, (lo, hi) in bounds.items():
        mask &= df[column].between(lo, hi)
    return df[mask].copy()


def summarize(df: pd.DataFrame, bounds: dict):
    """How much each filter removes on its own, before combining them."""
    for column, (lo, hi) in bounds.items():
        kept = df[column].between(lo, hi).sum()
        print(f"{column} in [{lo}, {hi}]: keeps {kept}/{len(df)}")
    combined = apply_filters(df, bounds)
    print(f"Combined: keeps {len(combined)}/{len(df)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--bound", nargs=3, action="append", metavar=("COLUMN", "MIN", "MAX"),
        help="Override/add a filter bound, e.g. --bound rg_A 10 25. Repeatable.",
    )
    ap.add_argument("--dry-run", action="store_true", help="Only print summarize(), don't write passing.csv")
    ap.add_argument("--list-bounds", action="store_true", help="List filterable columns and their min/max, then exit")
    args = ap.parse_args()

    df = load_scores()

    if args.list_bounds:
        skip_cols = {"design_path", "config_id", "ss_A", "ss_B"}
        numeric_cols = [c for c in df.columns if c not in skip_cols and pd.api.types.is_numeric_dtype(df[c])]
        print("Filterable columns (column: min - max, n non-null):")
        for col in numeric_cols:
            data = df[col].dropna()
            if data.empty:
                print(f"  {col}: no data")
            else:
                print(f"  {col}: {data.min():.4g} - {data.max():.4g}  (n={len(data)})")
        return

    bounds = dict(BOUNDS)
    if args.bound:
        for column, lo, hi in args.bound:
            bounds[column] = (float(lo), float(hi))

    summarize(df, bounds)

    if args.dry_run:
        print("\033[93m(dry run: passing.csv not written)\033[0m")
        return

    passing = apply_filters(df, bounds)
    out_path = SCRIPT_DIR / "passing.csv"
    passing.to_csv(out_path, index=False)
    print(f"Wrote {len(passing)} passing designs to {out_path}")


if __name__ == "__main__":
    main()

