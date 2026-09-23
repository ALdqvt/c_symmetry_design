import sys

sys.path.insert(0, "/home/panda/Resources/software/rfd2")
import pickle
from pathlib import Path
import numpy as np

_SEQ_KEYS = ("contigmap_hal", "contigmap_ref", "inpaint_seq", "inpaint_str")
_SEQ_REQUIRED = frozenset(_SEQ_KEYS)


def get_trb(directory, config_id, design_n) -> dict:
    """
    Returns a trb_data object
    """
    # Combine the directory and config_id to find .trb file
    trb_path = Path(directory) / config_id / f"{config_id}_design_{design_n}.trb"
    # Unpickle .trb file
    with open(trb_path, "rb") as f:
        trb_data = pickle.load(f)
    # Return trb object
    return trb_data


def get_sequence_arrays(trb: dict):
    """
    Input: trb_data from unpickled .trb file.
    Returns: Four arrays with residue-informaiton, in this order:
        contigmap_hal, contigmap_ref, inpaint_seq, inpaint_str

    Raises KeyError if any of these keys is missing from the .trb data.
    """

    if not trb.keys() >= _SEQ_REQUIRED:
        raise KeyError(
            f"Missing keys in trb:{trb}\nKeys missing:{_SEQ_REQUIRED - trb.keys()}"
        )

    return tuple(np.asarray(trb[k]) for k in _SEQ_KEYS)


def inpaint_lookup(trb):
    """
    Maps (chain, resnum) -> bool using inpaint_seq convention:
    True = sequence identity KEPT
    False = masked/inpainted

    Assumption: #TODO: Verify
    contigmap_hal is longer than inpaint_seq because it includes the
    ligand's atomized entries at the tail. Truncating to inpaint_seq's
    length keeps all protein positions and drops only the ligand tail.
    """
    contigmap_hal, _, inpaint_seq, _ = get_sequence_arrays(trb)

    min_len = min(len(contigmap_hal), len(inpaint_seq))
    lookup = {}
    for i in range(min_len):
        key = (str(contigmap_hal[i][0]), int(contigmap_hal[i][1]))
        lookup[key] = bool(inpaint_seq[i])

    return lookup
