# main script: config -> place -> symmetrize -> write + manifest

from pathlib import Path
from Bio.PDB import PDBParser, PDBIO
import csv
from datetime import datetime

from motif_config import MotifConfig, generate_configs, motif_extent_worst_case
from align_inputs import get_com
from placement import place_motif
from symmetrise import symmetrise

ALIGNED_DIR = Path("input_structures/aligned")
OUTPUT_DIR = Path("input_structures/symmetrised_test")
OUTPUT_DIR.mkdir(exist_ok=True)

parser = PDBParser(QUIET=True)
io = PDBIO()


aligned_files = sorted(ALIGNED_DIR.glob("*_aligned.pdb"))
structure = parser.get_structure(aligned_files[0].stem, aligned_files[0])
chain_atoms = list(structure[0]["A"].get_atoms())


motif_radius = motif_extent_worst_case(chain_atoms)
r_min = 2 * motif_radius
r_max = r_min + motif_radius  # might change


manifest_path = OUTPUT_DIR / "manifest.csv"
manifest_exists = manifest_path.exists()

manifest_file = open(manifest_path, "a", newline="")
writer = csv.writer(manifest_file)
if not manifest_exists:
    writer.writerow(["filepath", "r", "theta", "phi", "zeta", "timestamp"])

configs = list(generate_configs(r_min=r_min, r_max=r_max, r_step=2.0, grid_step_deg=15))
print(f"About to generate {len(configs)} structures.")

failed = []

for config in generate_configs(r_min=r_min, r_max=r_max, r_step=2.0, grid_step_deg=15):
    try:
        placed = place_motif(structure, config)
        symmetrised = symmetrise(placed, n_copies=6, include_ligand=True)

        out_path = OUTPUT_DIR / (
            f"{aligned_files[0].stem}_r{config.r:.1f}_theta{config.theta:.1f}"
            f"_phi{config.phi:.1f}_zeta{config.zeta:.1f}.pdb"
        )
        io.set_structure(symmetrised)
        io.save(str(out_path))

        writer.writerow([str(out_path), config.r, config.theta, config.phi, config.zeta, datetime.now().isoformat()])
        print(f"Wrote {out_path}")

    except Exception as e:
        print(f"FAILED: config={config} - {e}")
        failed.append(config)

manifest_file.close()
print(f"\nDone. {len(failed)} failures out of {len(configs)}.")
