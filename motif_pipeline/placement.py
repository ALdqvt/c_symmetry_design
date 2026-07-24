# rotate/translate motif per config

import numpy as np
from scipy.spatial.transform import Rotation
from Bio.PDB import PDBParser, PDBIO
import copy


def spherical_to_cartesian(r, theta_deg, phi_deg):
    """theta = longitude (0-360), phi = colatitude (0-180), both degrees."""
    theta = np.deg2rad(theta_deg)
    phi = np.deg2rad(phi_deg)
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return np.array([x, y, z])


def place_motif(structure, config):
    placed = copy.deepcopy(structure)
    atoms = list(placed.get_atoms())

    x_axis = np.array([1.0, 0.0, 0.0])

    # Step 1: rotate u (currently +x) to point in the theta/phi direction
    target_dir = spherical_to_cartesian(1.0, config.theta, config.phi)  # unit vector
    rot_direction = Rotation.align_vectors([target_dir], [x_axis])[0]

    for atom in atoms:
        atom.set_coord(rot_direction.apply(atom.get_coord()))

    # Step 2: spin by zeta about the new u direction
    rot_zeta = Rotation.from_rotvec(target_dir * np.deg2rad(config.zeta))

    for atom in atoms:
        atom.set_coord(rot_zeta.apply(atom.get_coord()))

    # Step 3: translate COM to (-r, 0, 0) — always along the fixed x-axis,
    # independent of theta/phi orientation
    translation = -x_axis * config.r
    for atom in atoms:
        atom.set_coord(atom.get_coord() + translation)

    return placed


########
if __name__ == "__main__":
    from pathlib import Path
    from motif_pipeline.motif_config import MotifConfig

    parser = PDBParser(QUIET=True)
    io = PDBIO()

    ALIGNED_DIR = Path("../input_structures/aligned")
    OUTPUT_DIR = Path("../input_structures/placed_test")
    OUTPUT_DIR.mkdir(exist_ok=True)

    aligned_files = sorted(ALIGNED_DIR.glob("*_aligned.pdb"))
    structure = parser.get_structure(aligned_files[0].stem, aligned_files[0])

# sanity check: theta=0, phi=90 -> along original +x, zeta=0 -> no spin
    config = MotifConfig(r=15.0, theta=0.0, phi=90.0, zeta=180.0)

    placed = place_motif(structure, config)

    out_path = OUTPUT_DIR / f"{aligned_files[0].stem}_placed.pdb"
    io.set_structure(placed)
    io.save(str(out_path))
