import yaml

from pathlib import Path
from copy import deepcopy

def load_base_config(path:Path) -> dict:
    with open(path) as file:
        return yaml.safe_load(file)


def set_nested(config: dict, dotted_key: str, value):
    """Set config['a']['b']['c'] via 'a.b.c', creating missing dicts along the way."""
    keys = dotted_key.split(".")
    d = config
    for k in keys[:-1]:
        if k not in d or not isinstance(d[k], dict):
            d[k] = {}
        d = d[k]
    d[keys[-1]] = value

def build_config(base_config: dict, overrides: dict) -> dict:
    config = deepcopy(base_config)
    for dotted_key, value in overrides.items():
        set_nested(config, dotted_key, value)
    return config


def write_config(config: dict, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False, default_flow_style=False)


if __name__ == "__main__":
    base = load_base_config(Path("configs/base.yaml"))

    overrides = {
        "inference.input_pdb": "/projects/arvid/c_symmetry_design/rfd2/in/test_output_1.pdb",
        "inference.num_designs": 2000,
        "inference.output_prefix": "/projects/arvid/c_symmetry_design/rfd2/out/tests/test_from_a",
        "contigmap.contigs": ["50,A1-209_B1-209"],  # placeholder, adjust to real A/B ranges
        "contigmap.has_termini": [True, True],
        "transforms.configs.RejectOutOfMemoryHazards.max_size": 100000, # match aa_ppi.yaml's inference default
        "contigmap.inpaint_seq": ["A1-43,B170-209"]
        # "ppi.hotspot_res": None,  # fill in once you decide which B residues matter
    }

    new_config = build_config(base, overrides)
    write_config(new_config, Path("configs/design_ab.yaml"))
    print("Wrote configs/design_ab.yaml")