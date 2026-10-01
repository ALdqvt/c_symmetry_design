"""Combine the two chains of an MPNN input into a single monomer sequence.

When diffusing an interface for a C-symmetry homo-oligomer, RFdiffusion output
contains two chains (A and B) that describe the same motif from opposite sides.
This module collapses that two-chain sequence into one sequence with a single
copy of the motif, using the .trb data to tell the residue classes apart.

Residue classes
---------------
Each residue is classified from the trb flags (inpaint_seq, inpaint_str):

    seq  (False, False)  the interface residues of each chain, kept as they are
    str  (False, True)   motif residues whose sequence identity was hidden
                         from the model. They are re-predicted, since they
                         contact the interface residues (seq) of the other
                         chain.
    m    (True,  True)   motif residues with fixed identity and structure

(True, False) is not expected and raises NotImplementedError.

Layout
------
The motif is present in both chains, split into an m part and a str part.
Which part comes first depends on the design, giving two layouts. In the
joined sequence the blocks appear as either:

    layout 1:  a_m   | a_str | a_seq | b_seq | b_str | b_m
    layout 2:  a_str | a_m   | a_seq | b_seq | b_m   | b_str

In both, a_seq and b_seq sit between the two motif halves, and the m blocks
of the two chains overlap, offset by the length of the str block. The motif
is rebuilt from that overlap, with str residues taking priority over any m
residue at the same position. This is the same for both layouts.

The output is a_seq + M + b_seq: the two interface segments (a_seq from
chain A, b_seq from chain B), flanking a single copy of the motif M. The
motif appears in both chains, so the two copies are merged into one, and its
position in the original joined sequence is not preserved.

Input requirements
------------------
- trb_data: the dict loaded from the .trb file (see mpnn.input_utils.get_trb)
- sequence: the MPNN sequence with the chain separator ':' already removed,
  so that its indices line up with the trb arrays

Usage
-----
From Python:

    from mpnn.input_utils import get_trb
    from mpnn.combine_chains import build_combined

    trb = get_trb("rfd2/out/", "0a60e948", 6)
    joined = fasta_sequence.replace(":", "")
    monomer = build_combined(trb, joined)

To inspect the intermediate pieces:

    from mpnn.combine_chains import combine_chains, rebuild_motif

    a_seq, a_str, a_m, b_seq, b_str, b_m = combine_chains(trb, joined)
    M = rebuild_motif(a_str, a_m, b_str, b_m)

Raises ValueError if the sequence still contains ':', is shorter than the trb
data, contains a chain other than A or B, or if the m blocks of the two chains
do not overlap consistently.
"""

from mpnn.input_utils import get_sequence_arrays



def _residues_to_str(res):
    """Join (resnum, i, aa) tuples into a string, ordered by i. """
    return "".join(aa for _, _, aa in sorted(res, key=lambda t: t[1]))

def combine_chains(trb_data: dict, sequence: str) -> tuple:
    if ":" in sequence:
        raise ValueError("sequence still contains ':'; remove the chain separator first")

    contigmap_hal, _, inpaint_seq, inpaint_str = get_sequence_arrays(trb_data)

    n = min(len(contigmap_hal), len(inpaint_seq), len(inpaint_str))
    if len(sequence) < n:
        raise ValueError(f"sequence has {len(sequence)} residues but trb data covers {n}")


    a_seq, a_str, a_m, b_seq, b_str, b_m = [], [], [], [], [], []

    # (seq_flag, str_flag) -> (chain A list, chain B list)
    buckets = {
        (True, True): (a_m, b_m),
        (False, True): (a_str, b_str),
        (False, False): (a_seq, b_seq),
    }

    # zip truncates to shortest list, this is by design
    for i, (tpl, seq_flag, str_flag) in enumerate(
        zip(contigmap_hal, inpaint_seq, inpaint_str)
    ):
        chain, resnum = tpl
        if chain not in ('A', 'B'):
            raise ValueError(f"Unknown chain {chain!r} at index {i}")

        key = (seq_flag, str_flag)
        if key not in buckets:
            raise NotImplementedError(
                f"Unexpected flag combination at index {i}: "
                f"seq_flag={seq_flag}, str_flag={str_flag}"
            )

        a_list, b_list = buckets[key]
        (a_list if chain == 'A' else b_list).append((resnum, i, sequence[i]))

    return a_seq, a_str, a_m, b_seq, b_str, b_m

def rebuild_motif(a_str, a_m, b_str, b_m):
    """Rebuild motif M from the chain A and chain B pieces.

    Each argument is a list of (resnum, i, aa).
    The m blocks are aligned by finding where one starts inside the other,
    then the str blocks are placed next to them and overwrite any m residue.
    Works for both layouts (str before or after m in each chain).
    """
    if not a_m or not b_m:
        raise ValueError("Both chains need at least one m residue to align the motif")

    a_ms, b_ms = _residues_to_str(a_m), _residues_to_str(b_m)

    # Offset of each m block within M
    k = a_ms.find(b_ms[:20])
    if k >= 0:
        a_off, b_off = 0, k
    else:
        k = b_ms.find(a_ms[:20])
        if k < 0:
            raise ValueError("m blocks of chain A and B do not overlap")
        a_off, b_off = k, 0

    motif = {}  # position in M -> amino acid

    # m residues: the overlapping region must agree between chains
    for off, ms in ((a_off, a_ms), (b_off, b_ms)):
        for j, aa in enumerate(ms):
            pos = off + j
            if pos in motif and motif[pos] != aa:
                raise ValueError(f"m residues disagree at motif position {pos}")
            motif[pos] = aa

    # str residues: positioned relative to their own chain's m block, and they win
    str_positions = set()

    def place_str(str_res, m_res, off):
        m_start = min(i for _, i, _ in m_res)
        for _, i, aa in str_res:
            pos = off + i - m_start
            if pos in str_positions:
                raise ValueError(f"str residues from A and B collide at motif position {pos}")
            str_positions.add(pos)
            motif[pos] = aa

    place_str(a_str, a_m, a_off)
    place_str(b_str, b_m, b_off)

    positions = sorted(motif)
    if positions != list(range(positions[0], positions[0] + len(positions))):
        raise ValueError("Rebuilt motif has gaps")

    return "".join(motif[p] for p in positions)


def combine(a_seq, M, b_seq):
    return _residues_to_str(a_seq) + M + _residues_to_str(b_seq)


def build_combined(trb_data: dict, sequence: str) -> str:
    a_seq, a_str, a_m, b_seq, b_str, b_m = combine_chains(trb_data, sequence)
    M = rebuild_motif(a_str, a_m, b_str, b_m)
    return combine(a_seq, M, b_seq)






