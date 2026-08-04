# build_rfd2_configs.py
# Reads a motif_pipeline manifest.csv, filters to viable configs, and writes
# one RFD2 YAML config per config_id into configs/generated/.

import pandas as pd
from pathlib import Path

from rfd2.rfd2_config import load_base_config, build_config, write_config


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

BASE_CONFIG_PATH = SCRIPT_DIR / "configs" / "base.yaml"
GENERATED_CONFIG_DIR = SCRIPT_DIR / "configs" / "generated"
OUTPUT_DIR = PROJECT_DIR / "rfd2" / "out"

# TODO: generalize to loop over all output_structures/*/manifest.csv once
# multiple structures are in play. Hardcoded to one structure for now.
MANIFEST_PATH = PROJECT_DIR / "output_structures" / "silicatein_second_batch_01_aligned" / "manifest.csv"


FLANK_INTERFACE = 50   # residues built at the "close" termini, to bridge the interface
INPAINT_WINDOW = 40 # TODO
N_DESIGNS = 10

# Sanity filter: orientation_dist shouldn't wildly exceed what FLANK_INTERFACE
# residues can plausibly bridge. This is a rough guess -- inspect the real
# distribution of orientation_dist in your manifest (manifest.orientation_dist.describe())
# before trusting this cutoff for a real run.
MAX_BRIDGE_DIST = 30.0  # Angstroms

# Set to an integer to only build configs for a quick test run;
# set to None to process everything that passes filtering.
LIMIT = None 

##### ##### #####

def get_ids_in_window(chain_start, chain_end, terminus, window):
    """Returns the residue number range within `window` of the given
    terminus ('N' or 'C'), clipped to the chain's actual bounds."""
    if terminus == "N":
        lo, hi = chain_start, min(chain_start + window - 1, chain_end)
    else:  # "C"
        lo, hi = max(chain_end - window + 1, chain_start), chain_end
    return lo, hi


def build_inpaint_seq(chain_start, chain_end, orientation, window=INPAINT_WINDOW):
    """
    Masks sequence identity for ALL residues within `window` of the
    interface-facing terminus on each chain -- matching whichever side
    is being extended by diffusion (same orientation convention as build_contig).
    """
    if orientation == "A_N-B_C":
        a_lo, a_hi = get_ids_in_window(chain_start, chain_end, "N", window)
        b_lo, b_hi = get_ids_in_window(chain_start, chain_end, "C", window)
    elif orientation == "A_C-B_N":
        a_lo, a_hi = get_ids_in_window(chain_start, chain_end, "C", window)
        b_lo, b_hi = get_ids_in_window(chain_start, chain_end, "N", window)
    else:
        raise ValueError(f"Unexpected orientation value: {orientation!r}")

    a_range = f"A{a_lo}-{a_hi}" if a_lo != a_hi else f"A{a_lo}"
    b_range = f"B{b_lo}-{b_hi}" if b_lo != b_hi else f"B{b_lo}"
    return f"{a_range},{b_range}"

def build_contig(chain_start, chain_end, orientation):
    """
    Builds a contig string with asymmetric flanks: a large flank at the
    termini that are close together (per orientation), and a small flank
    at the far termini.

    orientation: 'A_N-B_C' or 'A_C-B_N', from decide_interface_orientation.
        'A_N-B_C' -> A's N-terminus is close to B's C-terminus.
        'A_C-B_N' -> A's C-terminus is close to B's N-terminus.

    Contig segment format is {N-side flank},{chain}{start}-{end},{C-side flank}
    i.e. the first number is always the N-terminal side, the second is
    always the C-terminal side, regardless of which side is "interface".
    """
    if orientation not in ("A_N-B_C", "A_C-B_N"):
        raise ValueError(f"Unexpected orientation value: {orientation!r}")

    if orientation == "A_N-B_C":
        # A's N-terminus and B's C-terminus are the close pair -> grow there
        a_seg = f"{FLANK_INTERFACE},A{chain_start}-{chain_end}"
        b_seg = f"B{chain_start}-{chain_end},{FLANK_INTERFACE}"
    else:  # "A_C-B_N"
        a_seg = f"A{chain_start}-{chain_end},{FLANK_INTERFACE}"
        b_seg = f"{FLANK_INTERFACE},B{chain_start}-{chain_end}"

    return f"{a_seg}_{b_seg}"


def load_and_filter_manifest(manifest_path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path)

    n_total = len(manifest)
    n_success = (manifest["status"] == "success").sum()

    filtered = manifest[manifest["status"] == "success"].copy()
    filtered = filtered[filtered["orientation_dist"] < MAX_BRIDGE_DIST]

    print(
        f"Manifest: {n_total} total rows, {n_success} successful, "
        f"{len(filtered)} pass orientation_dist < {MAX_BRIDGE_DIST}."
    )

    if LIMIT is not None:
        filtered = filtered.head(LIMIT)
        print(f"LIMIT set -- only processing first {len(filtered)} rows.")

    return filtered


def main():
    GENERATED_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    base = load_base_config(BASE_CONFIG_PATH)

    rows = load_and_filter_manifest(MANIFEST_PATH)

    n_written = 0
    n_failed = 0

    for _, row in rows.iterrows():
        config_id = row["config_id"]
        try:
            contig = build_contig(row["chain_start"], row["chain_end"], row["orientation"])
            inpaint_seq_value = build_inpaint_seq(row["chain_start"], row["chain_end"], row["orientation"])

            overrides = {
                "inference.input_pdb": str((PROJECT_DIR / row["filepath"]).resolve()),
                "inference.output_prefix": str(OUTPUT_DIR / config_id / f"{config_id}_design"),
                "inference.num_designs": N_DESIGNS,
                "inference.ligand": "LIG",
                "contigmap.contigs": [contig],
                "contigmap.has_termini": [True, True],
                "contigmap.inpaint_seq": [inpaint_seq_value],
                "transforms.configs.RejectOutOfMemoryHazards.max_size": 100000,
            }

            new_config = build_config(base, overrides)
            write_config(new_config, GENERATED_CONFIG_DIR / f"{config_id}.yaml")
            n_written += 1

        except Exception as e:
            print(f"FAILED to build config for config_id={config_id}: {e}")
            n_failed += 1

    print(f"\nDone. Wrote {n_written} configs to {GENERATED_CONFIG_DIR}, {n_failed} failed.")


if __name__ == "__main__":
    main()
