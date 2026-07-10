import numpy as np
from scipy.spatial.transform import Rotation
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO




def get_com(atoms):
    coords = np.array([a.get_coord() for a in atoms])
    return coords.mean(axis=0)

def align_structure(structure, n_term_atom, c_term_atom):
    atoms = list(structure.get_atoms())
    a_chain = structure[0]["A"]
    b_chain = structure[0]["B"]
    a_chain_atoms = list(a_chain.get_atoms())
    com = get_com(a_chain_atoms)


    # Step 1: translate COM to origin
    for atom in atoms:
        atom.set_coord(atom.get_coord() - com)

    ligand_atoms = list(b_chain.get_atoms())
    lig_coord = get_com(ligand_atoms)

    # Step 2: rotate u -> x-axis
    u = lig_coord / np.linalg.norm(lig_coord)
    x_axis = np.array([1.0, 0.0, 0.0])
    rot1 = Rotation.align_vectors([x_axis], [u])[0]

    for atom in atoms:
        atom.set_coord(rot1.apply(atom.get_coord()))

    # Step 3: rotate about x-axis so N->C projects onto +y as closely as possible
    n_coord = n_term_atom.get_coord()
    c_coord = c_term_atom.get_coord()
    nc_vec = c_coord - n_coord
    nc_proj = nc_vec.copy()
    nc_proj[0] = 0  # drop x-component, project onto y-z plane

    if np.linalg.norm(nc_proj) < 1e-6:
        # N->C is (near) parallel to u -- projection undefined, skip step 3
        rot2 = Rotation.identity()
    else:
        nc_proj /= np.linalg.norm(nc_proj)
        y_axis = np.array([0.0, 1.0, 0.0])
        angle = np.arctan2(
            np.dot(np.cross(nc_proj, y_axis), x_axis),
            np.dot(nc_proj, y_axis)
        )
        rot2 = Rotation.from_rotvec(x_axis * angle)

    for atom in atoms:
        atom.set_coord(rot2.apply(atom.get_coord()))

    return structure

def main():
    INPUT_DIR = Path("input_structures")
    OUTPUT_DIR = Path("input_structures/aligned")
    OUTPUT_DIR.mkdir(exist_ok=True)


    parser = PDBParser(QUIET=True)
    io = PDBIO()

    failed = []
    pdb_files = sorted(INPUT_DIR.glob("*.pdb"))

    for pdb_path in pdb_files:
        try:
            structure = parser.get_structure(pdb_path.stem, pdb_path)
            a_chain = structure[0]["A"]

            residue_ids = sorted(res.id[1] for res in a_chain if res.id[0] == " ")
            n_term_atom = a_chain[residue_ids[0]]["CA"]
            c_term_atom = a_chain[residue_ids[-1]]["CA"]

            aligned = align_structure(structure, n_term_atom, c_term_atom)

            out_path = OUTPUT_DIR / f"{pdb_path.stem}_aligned.pdb"
            io.set_structure(aligned)
            io.save(str(out_path))

        except Exception as e:
            print(f"FAILED: {pdb_path.name} — {e}")
            failed.append(pdb_path.name)

    print(f"\nDone. {len(failed)} failures out of {len(pdb_files)}")
    if failed:
        print("Failed files:", failed)

if __name__ == "__main__":
    main()