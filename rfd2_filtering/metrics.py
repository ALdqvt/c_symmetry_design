"""Structural quality metrics for filtering RFD2 output designs.

Two metrics implemented now, both pure geometry so they need no external
tools:
  - clash score: non-bonded atom pairs closer than threshold, within and
    across chains
  - backbone continuity: flags chain breaks / implausibly stretched
    backbone via consecutive CA-CA distance

Stubs below (to fill in later):
  - secondary_structure: DSSP-based helix/sheet/loop fractions
  - intertwining: topological interlocking check between chains
"""

import numpy as np
from scipy.spatial import cKDTree
import warnings # Suppress warnings about mmCIF format from DSSP as we're dealing with PDB files
warnings.filterwarnings("ignore", message="parse error at line 1", module="Bio.PDB.DSSP")

# Typical peptide bond CA-CA distance is ~3.8 A. Flag anything notably
# longer as a likely chain break, or an unrealistically stretched-out
# backbone -- this also catches the "too-long residue distances" you
# noticed, since that shows up as an inflated gap here too. Tune once
# you've looked at the real distribution (see filtering.plot_metric).
BACKBONE_BREAK_THRESHOLD = 4.0  # Angstroms, CA(i) -> CA(i+1)
# Yet to find a good balance on bb_break_thrshld. Anything above 4.0 Å is not normal though so
# it's a valid threshold for removing chain breaks.

CLASH_THRESHOLD = 2.0  # Angstroms, matches interface_checks.check_clash default
CLASH_MIN_SEQ_SEP = 2  # ignore pairs this close in sequence on the same chain (bonded/near-bonded)

# mkdssp isn't installable into ca_rfd (Python 3.9, no compatible conda-forge
# build) so we call the binary from the c-symmetry-design env directly,
# regardless of which env this script is actually running under.
MKDSSP_PATH = "/home/panda/anaconda3/envs/c-symmetry-design/bin/mkdssp"

def get_all_atoms(structure, exclude_hetero=True):
    atoms = []
    for atom in structure.get_atoms():
        residue = atom.get_parent()
        if exclude_hetero and residue.id[0] != " ":
            continue
        atoms.append(atom)
    return atoms


def compute_clash_score(structure, threshold=CLASH_THRESHOLD, min_seq_sep=CLASH_MIN_SEQ_SEP):
    """
    Counts non-bonded atom pairs closer than `threshold` Angstroms, within
    and across chains. Pairs within `min_seq_sep` residues of each other on
    the same chain are skipped, since those are expected to be close.

    Returns {"n_clashes": int, "min_clash_dist": float | None}.
    """
    atoms = get_all_atoms(structure)
    if len(atoms) < 2:
        return {"n_clashes": 0, "min_clash_dist": None}

    coords = np.array([a.get_coord() for a in atoms])
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=threshold, output_type="ndarray")

    n_clashes = 0
    min_dist = None
    for i, j in pairs:
        res_i = atoms[i].get_parent()
        res_j = atoms[j].get_parent()
        same_chain = res_i.get_parent().id == res_j.get_parent().id
        if same_chain and abs(res_i.id[1] - res_j.id[1]) < min_seq_sep:
            continue
        dist = float(np.linalg.norm(coords[i] - coords[j]))
        n_clashes += 1
        if min_dist is None or dist < min_dist:
            min_dist = dist

    return {"n_clashes": n_clashes, "min_clash_dist": min_dist}


def compute_backbone_continuity(structure, break_threshold=BACKBONE_BREAK_THRESHOLD):
    """
    Walks each chain's CA trace in residue order and flags gaps bigger than
    `break_threshold`.

    Returns {"max_ca_gap": float, "n_backbone_breaks": int}.
    """
    max_gap = 0.0
    n_breaks = 0

    for chain in structure.get_chains():
        residues = [r for r in chain if r.id[0] == " " and "CA" in r]
        residues.sort(key=lambda r: r.id[1])
        for r1, r2 in zip(residues, residues[1:]):
            gap = float(np.linalg.norm(r1["CA"].get_coord() - r2["CA"].get_coord()))
            max_gap = max(max_gap, gap)
            if gap > break_threshold:
                n_breaks += 1

    return {"max_ca_gap": round(max_gap, 3), "n_backbone_breaks": n_breaks}


def compute_secondary_structure(pdb_path, trb_data):
    """
    Runs DSSP (via Biopython's wrapper around the `mkdssp` binary) and
    returns:
      - the full per-residue secondary structure string per chain, for
        closer inspection later
      - helix/sheet/loop fractions computed ONLY over the diffused
        extension residues (i.e. excluding the original motif), since the
        motif's structure is already fixed/known and would otherwise
        dilute the signal on what diffusion actually built

    DSSP per-residue codes: H/G/I -> helix (h), E/B -> sheet (e), anything
    else -> loop/coil/turn (l). The ss string uses these single-letter
    codes in residue order, one string per chain.

    Only meant to be run on designs that already passed clash and
    backbone-continuity filtering -- DSSP is comparatively slow, and a
    clashing or broken backbone isn't worth spending it on.

    Returns:
        {
            "ss_A": "lllhhhhhhh...", "ss_B": "...",
            "frac_helix": float, "frac_sheet": float, "frac_loop": float,
        }
    frac_* are computed over extension residues across all chains combined
    (not per-chain), matching the existing HEADER columns.
    """
    import tempfile
    from Bio.PDB import PDBParser, PDBIO, DSSP

    structure = PDBParser(QUIET=True).get_structure(pdb_path.stem, pdb_path)

    motif_by_chain = {}
    for chain_id, resnum in trb_data["con_hal_pdb_idx"]:
        motif_by_chain.setdefault(chain_id, set()).add(resnum)

    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w") as tmp:
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp.name)

        model = structure[0]
        dssp = DSSP(model, tmp.name, dssp = MKDSSP_PATH)

        ss_by_chain = {}
        n_helix = n_sheet = n_loop = 0  # extension-only counts

        for key in dssp.keys():
            chain_id, res_id = key[0], key[1]
            resnum = res_id[1]
            ss_code = dssp[key][2]

            if ss_code in ("H", "G", "I"):
                letter, bucket = "h", "helix"
            elif ss_code in ("E", "B"):
                letter, bucket = "e", "sheet"
            else:
                letter, bucket = "l", "loop"

            ss_by_chain.setdefault(chain_id, []).append(letter)

            is_motif = resnum in motif_by_chain.get(chain_id, set())
            if not is_motif:
                if bucket == "helix":
                    n_helix += 1
                elif bucket == "sheet":
                    n_sheet += 1
                else:
                    n_loop += 1

    result = {f"ss_{cid}": "".join(letters) for cid, letters in ss_by_chain.items()}

    n_total = n_helix + n_sheet + n_loop
    if n_total == 0:
        result.update({"frac_helix": None, "frac_sheet": None, "frac_loop": None})
    else:
        result.update({
            "frac_helix": round(n_helix / n_total, 4),
            "frac_sheet": round(n_sheet / n_total, 4),
            "frac_loop": round(n_loop / n_total, 4),
        })

    return result

def compute_radius_of_gyration(structure):
    """
    Per-chain radius of gyration (CA-only). No normalization -- look at
    the distribution across your batch (same as the other metrics) and
    flag outliers directly, rather than comparing to a theoretical
    expected value.

    Returns {"rg_A": ..., "rg_B": ...} for whichever chains are present.
    """
    result = {}
    for chain in structure.get_chains():
        residues = [r for r in chain if r.id[0] == " " and "CA" in r]
        if len(residues) < 3:
            continue
        coords = np.array([r["CA"].get_coord() for r in residues])
        com = coords.mean(axis=0)
        rg = float(np.sqrt(np.mean(np.sum((coords - com) ** 2, axis=1))))
        result[f"rg_{chain.id}"] = round(rg, 3)
    return result

def compute_com_shift(structure, trb_data):
    """
    For each chain, compares the center of mass computed over only the
    original motif residues (non-diffused, per RFD2's own bookkeeping)
    against the center of mass computed over the full chain (motif +
    diffused extension). A large shift means the extension is dragging
    that chain's mass distribution away from where the original motif
    sat -- e.g. an interface extension reaching too far toward the other
    chain.

    trb_data: the dict loaded from the design's .trb file. Uses
    con_hal_pdb_idx, which lists (chain_id, resnum) pairs in the *output*
    (hallucinated) numbering for residues that came from the motif.

    Returns {"com_shift_A": ..., "com_shift_B": ...} for whichever chains
    have motif residues present.
    """
    motif_by_chain = {}
    for chain_id, resnum in trb_data["con_hal_pdb_idx"]:
        motif_by_chain.setdefault(chain_id, set()).add(resnum)

    result = {}
    for chain in structure.get_chains():
        cid = chain.id
        if cid not in motif_by_chain:
            continue
        motif_resnums = motif_by_chain[cid]

        all_residues = [r for r in chain if r.id[0] == " " and "CA" in r]
        motif_residues = [r for r in all_residues if r.id[1] in motif_resnums]
        if not motif_residues or not all_residues:
            continue

        motif_coords = np.array([r["CA"].get_coord() for r in motif_residues])
        full_coords = np.array([r["CA"].get_coord() for r in all_residues])

        shift = float(np.linalg.norm(motif_coords.mean(axis=0) - full_coords.mean(axis=0)))
        result[f"com_shift_{cid}"] = round(shift, 3)

    return result

# --- Stubs for future metrics ---
def compute_intertwining(structure):
    """
    TODO: flag good-interface-but-topologically-interlocked designs, i.e.
    chains that look like a clean interface but are mechanically locked
    together (can't be pulled apart without breaking a bond). A cheap first
    pass: check whether one chain's CA trace threads through loop closures
    formed by the other chain; a more rigorous version would compute a
    linking number between the two backbone curves.
    """
    raise NotImplementedError
