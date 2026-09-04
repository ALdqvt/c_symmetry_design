"""
Walks RFD2 output PDBs and writes a per-design scores CSV.

Usage:
    python -m rfd2_filtering.run_scoring \
        --out-dir rfd2/out --scores-path rfd2_filtering/scores.csv

See the NOTE at the bottom for hooking this into run_all.sh as an
"on the go" step once you're ready for that.
"""
import sys 
sys.path.insert(0, "/home/panda/Resources/software/rfd2") # Needed for unpickling .trb files

import pickle
import argparse
import csv
from pathlib import Path

from Bio.PDB import PDBParser

from rfd2_filtering.metrics import (
    compute_clash_score,
    compute_backbone_continuity,
    compute_radius_of_gyration,
    compute_com_shift,
    compute_secondary_structure,
)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

HEADER = [
    "design_path",
    "config_id",
    "n_clashes",
    "min_clash_dist",
    "max_ca_gap",
    "n_backbone_breaks",
    "frac_helix",
    "frac_sheet",
    "frac_loop",
    "ss_A",
    "ss_B",
    "rg_A",
    "rg_B",
    "com_shift_A",
    "com_shift_B",
    # intertwining columns get appended here once
    # compute_intertwining are implemented
]

_parser_pdb = PDBParser(QUIET=True)


def score_design(pdb_path: Path, config_id: str) -> dict:
    structure = _parser_pdb.get_structure(pdb_path.stem, pdb_path)
    row = {"design_path": str(pdb_path), "config_id": config_id}
    row.update(compute_clash_score(structure))
    row.update(compute_backbone_continuity(structure))
    row.update(compute_radius_of_gyration(structure))

    trb_path = pdb_path.with_suffix(".trb")
    if trb_path.exists():
        with open(trb_path, "rb") as f:
            trb_data = pickle.load(f)
        row.update(compute_com_shift(structure, trb_data))
        if row.get("n_clashes") == 0 and row.get("n_backbone_breaks") == 0:
            row.update(compute_secondary_structure(pdb_path, trb_data))
    else:
        print(f"WARNING: no .trb found for {pdb_path}, skipping com_shift")

    

    return row


def find_designs(out_dir: Path):
    """
    RFD2 writes designs as {output_prefix}_{i}.pdb, and output_prefix is
    OUTPUT_DIR/config_id/design (see build_rfd2_configs.py), so each
    config_id gets its own subdirectory under out_dir.

    NOTE: confirm "design_*.pdb" matches your actual RFD2 output naming --
    adjust the glob below if it writes a different suffix/extension.
    """
    for config_dir in sorted(p for p in out_dir.iterdir() if p.is_dir()):
        for pdb_path in sorted(config_dir.glob(f"{config_dir.name}_design_*.pdb")):
            yield config_dir.name, pdb_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=PROJECT_DIR / "rfd2" / "out")
    ap.add_argument("--scores-path", type=Path, default=SCRIPT_DIR / "scores.csv")
    args = ap.parse_args()

    args.scores_path.parent.mkdir(parents=True, exist_ok=True)

    n_scored, n_failed = 0, 0
    current_config = None
    with open(args.scores_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADER)
        writer.writeheader()
        for config_id, pdb_path in find_designs(args.out_dir):
            if config_id != current_config:
                print(f"Scoring {config_id} ...")
                current_config = config_id
            try:
                row = score_design(pdb_path, config_id)
                writer.writerow(row)
                n_scored += 1
            except Exception as e:
                print(f"FAILED to score {pdb_path}: {e}")
                n_failed += 1

    print(f"Scoring complete: {n_scored} designs scored, ({n_failed} unable to score)")
    print(f"Wrote scores to {args.scores_path}")


if __name__ == "__main__":
    main()

# NOTE on "on the go" filtering in run_all.sh:
# This only needs an --out-dir of freshly written PDBs, so it can run as a
# step straight after the RFD2 inference call, e.g. appending:
#   python -m rfd2_filtering.run_scoring --out-dir rfd2/out --scores-path rfd2_filtering/scores.csv
# Threshold filtering (filtering.py) is a separate, deliberately interactive
# step meant to be driven by a human looking at plots -- not auto-run in the
# pipeline. If you later want a fully automated pass/fail gate, we can add a
# --thresholds file to filtering.py and call it non-interactively too.
