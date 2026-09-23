import json
from pathlib import Path
from datetime import datetime

import pandas as pd

from mpnn.input_utils import get_trb, inpaint_lookup

SCRIPT_DIR = Path(__file__).resolve().parent


def build_redesigned_residues(trb) -> str:
    """
    Returns space separated string like "A7 A8 A9" of residues for redesign.
    """
    lookup = inpaint_lookup(trb)
    positions = [
        f"{chain}{resnum}"
        for (chain, resnum), fixed_residue in lookup.items()
        if not fixed_residue
    ]
    return " ".join(positions)


def main(pdb_dir, passing_csv_path):
    out_dir = SCRIPT_DIR / "input_files" / datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(passing_csv_path, dtype={"config_id": str, "design_path": str})

    pdb_paths = []
    redesigned = {}

    for _, row in df.iterrows():
        pdb_path = str(Path(str(row["design_path"])).resolve())
        if not Path(pdb_path).exists():
            raise FileNotFoundError(f"Resolved path doesn not exist: {pdb_path}")
        config_id = row["config_id"]
        design_n = int(Path(pdb_path).stem.rsplit("_design_", 1)[1])

        trb = get_trb(pdb_dir, config_id, design_n)

        pdb_paths.append(pdb_path)
        redesigned[pdb_path] = build_redesigned_residues(trb)

    with open(out_dir / "pdb_ids.json", "w") as f:
        json.dump(pdb_paths, f, indent=2)

    with open(out_dir / "redesigned_residues_multi.json", "w") as f:
        json.dump(redesigned, f, indent=2)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--passing-csv-path", required=True)
    args = ap.parse_args()

    main(args.pdb_dir, args.passing_csv_path)
