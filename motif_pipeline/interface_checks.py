import numpy as np
from scipy.spatial import cKDTree


def check_clash(atoms_a, atoms_b, threshold=2.0):

    coords_a = np.array([a.get_coord() for a in atoms_a])
    coords_b = np.array([a.get_coord() for a in atoms_b])
    tree_b = cKDTree(coords_b)
    dists, _ = tree_b.query(coords_a, k=1)
    min_dist = float(np.min(dists))
    return min_dist < threshold, min_dist


def get_terminus_ca(chain, terminus="N"):
    first, last = get_residue_range(chain)
    res_id = first if terminus == "N" else last
    return chain[res_id]["CA"].get_coord()


def decide_interface_orientation(chain_a, chain_b):

    a_n = get_terminus_ca(chain_a, "N")
    a_c = get_terminus_ca(chain_a, "C")
    b_n = get_terminus_ca(chain_b, "N")
    b_c = get_terminus_ca(chain_b, "C")

    dist_an_bc = np.linalg.norm(a_n - b_c)  # A's N-term near B's C-term
    dist_ac_bn = np.linalg.norm(a_c - b_n)  # A's C-term near B's N-term

    if dist_an_bc < dist_ac_bn:
        return "A_N-B_C", dist_an_bc, dist_ac_bn
    else:
        return "A_C-B_N", dist_ac_bn, dist_an_bc


def get_residue_range(chain):
    """Returns (first_res, last_resi) for a chain, ignoring heteroatoms"""
    residue_ids = sorted(res.id[1] for res in chain if res.id[0] == " ")
    return residue_ids[0], residue_ids[-1]
