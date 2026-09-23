"""
Grafts the full A-M-B monomer from each passing design (chain A's
extension superimposed onto chain B's frame), checks the graft for
clashes, and writes clean monomer PDBs ready for MPNN.

Usage:
    python -m rfd2_filtering.graft_monomers \
        --passing-path rfd2_filtering/passing.csv \
        --monomers-dir rfd2_filtering/monomers
"""

import sys

sys.path.insert(
    0, "/home/panda/Resources/software/rfd2"
)  # needed to unpickle .trb files

import argparse
import pickle
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser, PDBIO

from rfd2_filtering.metrics import build_grafted_monomer

SCRIPT_DIR = Path(__file__).resolve().parent

_parser_pdb = PDBParser(QUIET=True)
_io = PDBIO()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--passing-path", type=Path, default=SCRIPT_DIR / "passing.csv")
    ap.add_argument("--monomers-dir", type=Path, default=SCRIPT_DIR / "monomers")
    args = ap.parse_args()

    df = pd.read_csv(args.passing_path, dtype={"config_id": str, "design_path": str})
    args.monomers_dir.mkdir(parents=True, exist_ok=True)

    n_written, n_failed = 0, 0
    for _, row in df.iterrows():
        pdb_path = Path(row["design_path"])
        trb_path = pdb_path.with_suffix(".trb")
        stem = pdb_path.stem

        try:
            structure = _parser_pdb.get_structure(stem, pdb_path)
            with open(trb_path, "rb") as f:
                trb_data = pickle.load(f)

            monomer, _, _ = build_grafted_monomer(structure, trb_data)

            out_path = args.monomers_dir / f"{stem}_monomer.pdb"
            _io.set_structure(monomer)
            _io.save(str(out_path))
            n_written += 1

        except Exception as e:
            print(f"FAILED to graft {stem}: {e}")
            n_failed += 1

    print(f"Grafting complete: {n_written} monomers written, {n_failed} failed.")
    print(f"Wrote monomers to {args.monomers_dir}")


if __name__ == "__main__":
    main()
