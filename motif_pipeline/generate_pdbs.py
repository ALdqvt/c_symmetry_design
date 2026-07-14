# main script: config -> place -> symmetrize -> write + manifest

from pathlib import Path
from Bio.PDB import PDBParser, PDBIO
import csv
from datetime import datetime

from motif_config import generate_configs, motif_extent_worst_case
from placement import place_motif
from symmetrise import symmetrise

ALIGNED_DIR = Path("../input_structures/aligned")
OUTPUT_ROOT = Path("output_structures")
OUTPUT_ROOT.mkdir(exist_ok=True)

R_STEP = 2.0
GRID_STEP_DEG = 15

parser = PDBParser(QUIET=True)
io = PDBIO()

aligned_files = sorted(ALIGNED_DIR.glob("*_aligned.pdb"))

for aligned_path in aligned_files:
    out_dir = OUTPUT_ROOT / aligned_path.stem
    out_dir.mkdir(exist_ok=True)

    structure = parser.get_structure(aligned_path.stem, aligned_path)
    chain_atoms = list(structure[0]["A"].get_atoms())

    motif_radius = motif_extent_worst_case(chain_atoms)
    r_min = 2 * motif_radius
    r_max = r_min + motif_radius  # might change

    manifest_path = out_dir / "manifest.csv"
    manifest_exists = manifest_path.exists()

    manifest_file = open(manifest_path, "a", newline="")
    writer = csv.writer(manifest_file)
    if not manifest_exists:
        writer.writerow(["filepath", "r", "theta", "phi", "zeta", "timestamp"])

    configs = list(generate_configs(r_min=r_min, r_max=r_max, r_step=R_STEP, grid_step_deg=GRID_STEP_DEG))
    print(f"About to generate {len(configs)} structures.")

    failed = []

    for config in configs:
        try:
            placed = place_motif(structure, config)
            symmetrised = symmetrise(placed, n_copies=6, include_ligand=True)

            out_path = out_dir / (
                f"{aligned_path.stem}_r{config.r:.1f}_theta{config.theta:.1f}"
                f"_phi{config.phi:.1f}_zeta{config.zeta:.1f}.pdb"
            )
            io.set_structure(symmetrised)
            io.save(str(out_path))

            writer.writerow(
                [str(out_path), config.r, config.theta, config.phi, config.zeta, datetime.now().isoformat()])
            print(f"Wrote {out_path}")

        except Exception as e:
            print(f"FAILED: config={config} - {e}")
            failed.append(config)

    manifest_file.close()
    print(f"\nDone. {len(failed)} failures out of {len(configs)}.")
