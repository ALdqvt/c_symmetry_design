# C6 replication

import numpy as np
import copy
import string
from scipy.spatial.transform import Rotation
from Bio.PDB import Structure, Model


def symmetrise(placed_structure, n_copies=6, include_ligand=False, lig_chain_id="B"):
    z_axis = np.array([0.0, 0.0, 1.0])
    chain_ids = string.ascii_uppercase[:n_copies]  # A, B, C, D, E, F

    new_structure = Structure.Structure("symmetrized")
    new_model = Model.Model(0)
    new_structure.add(new_model)

    # assume the motif chain in placed_structure is "A"
    source_chain = placed_structure[0]["A"]
    source_ligand = None
    if include_ligand:
        source_ligand = placed_structure[0][lig_chain_id]
    lig_chain_ids = string.ascii_uppercase[-n_copies:]

    for i, chain_id in enumerate(chain_ids):
        angle = 360.0 / n_copies * i
        rot = Rotation.from_rotvec(z_axis * np.deg2rad(angle))

        new_chain = copy.deepcopy(source_chain)
        new_chain.id = chain_id

        for atom in new_chain.get_atoms():
            atom.set_coord(rot.apply(atom.get_coord()))

        new_model.add(new_chain)

        if include_ligand:
            new_ligand = copy.deepcopy(source_ligand)
            new_ligand.id = lig_chain_ids[i]
            for atom in new_ligand.get_atoms():
                atom.set_coord(rot.apply(atom.get_coord()))
            new_model.add(new_ligand)

    return new_structure
