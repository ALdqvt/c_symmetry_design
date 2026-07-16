# main script: config -> place -> symmetrize -> write + manifest

from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Structure
import csv
from datetime import datetime
import hashlib

from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from motif_config import generate_configs, motif_extent_worst_case
from motif_pipeline.align_inputs import get_com
from placement import place_motif
from symmetrise import symmetrise

ALIGNED_DIR = Path("../input_structures/aligned")
OUTPUT_ROOT = Path("output_structures")
OUTPUT_ROOT.mkdir(exist_ok=True)

R_STEP = 2.0
GRID_STEP_DEG = 15

N_COPIES = 6
INCLUDE_LIGAND = True

HEADER = [
        "config_id", "structure_name", "filepath",
        "r", "theta", "phi", "zeta",
        "n_copies", "include_ligand",
        "r_min", "r_max", "motif_radius",
        "status", "error", "timestamp",
    ]



def make_config_id(structure_name, config):
    key = f"{structure_name}_{config.r:.2f}_{config.theta:.2f}_{config.phi:.2f}_{config.zeta:.2f}"
    return hashlib.sha1(key.encode()).hexdigest()[:8]

def add_ori_atom_to_struct(struct: Structure, com:np.array):
    ori_atom = Atom(name="ORI",
                    coord=com,
                    bfactor=0.0,
                    occupancy=0.0,
                    altloc=" ",
                    fullname=" ORI",
                    serial_number=99999,
                    element="XE")
    ori_residue = Residue(id=("H_ORI", 1, " "), resname="ORI", segid=" ")
    ori_residue.add(ori_atom)
    model = struct[0]
    if "P" in model:
        chain_p = model["P"]
    else:
        chain_p = Chain("P")
        model.add(chain_p)
    chain_p.add(ori_residue)
    return struct

parser = PDBParser(QUIET=True)
io = PDBIO()

aligned_files = sorted(ALIGNED_DIR.glob("*_aligned.pdb"))

for aligned_path in aligned_files:
    out_dir = OUTPUT_ROOT / aligned_path.stem
    out_dir.mkdir(exist_ok=True)


    structure = parser.get_structure(aligned_path.stem, aligned_path)
    chain_atoms = list(structure[0]["A"].get_atoms())

    motif_radius = motif_extent_worst_case(chain_atoms)
    r_min = 1.5 * motif_radius
    r_max = r_min + motif_radius  # might change

    manifest_path = out_dir / "manifest.csv"
    manifest_exists = manifest_path.exists()
    manifest_file = open(manifest_path, "a", newline="")
    writer = csv.writer(manifest_file)
    if not manifest_exists:
        writer.writerow(HEADER)

    configs = list(generate_configs(r_min=r_min, r_max=r_max, r_step=R_STEP, grid_step_deg=GRID_STEP_DEG))
    print(f"About to generate {len(configs)} structures.")

    failed = []

    for config in configs:
        config_id = make_config_id(aligned_path.stem, config)
        try:
            placed = place_motif(structure, config)
            symmetrised = symmetrise(placed, n_copies=N_COPIES, include_ligand=INCLUDE_LIGAND)
            chain_a = symmetrised[0]["A"]
            chain_b = symmetrised[0]["B"]
            a_and_b = list(chain_a.get_atoms())
            a_and_b.extend(list(chain_b.get_atoms()))
            com = get_com(a_and_b)
            symmetrised = add_ori_atom_to_struct(symmetrised, com)

            out_path = out_dir / (
                f"{aligned_path.stem}_r{config.r:.1f}_theta{config.theta:.1f}"
                f"_phi{config.phi:.1f}_zeta{config.zeta:.1f}.pdb"
            )
            io.set_structure(symmetrised)
            io.save(str(out_path))

            writer.writerow([
                config_id, aligned_path.stem, str(out_path),
                config.r, config.theta, config.phi, config.zeta,
                N_COPIES, INCLUDE_LIGAND,
                r_min, r_max, motif_radius,
                "success", "", datetime.now().isoformat(),
            ])
            print(f"Wrote {out_path}")

        except Exception as e:
            writer.writerow([
                config_id, aligned_path.stem, "",
                config.r, config.theta, config.phi, config.zeta,
                N_COPIES, INCLUDE_LIGAND,
                r_min, r_max, motif_radius,
                "failed", str(e), datetime.now().isoformat(),
            ])
            print(f"FAILED: config={config} - {e}")
            failed.append(config)

    manifest_file.close()
    print(f"\nDone. {len(failed)} failures out of {len(configs)}.")
