# How to load a trb file:

```bash
conda activate ca_rfd
cd /projects/arvid/c_symmetry_design
export PYTHONPATH=/home/panda/Resources/software/rfd2
python
```


```python
import sys
sys.path.insert(0, "/home/panda/Resources/software/rfd2")

import pickle
from Bio.PDB import PDBParser
from rfd2_filtering.metrics import build_grafted_monomer

config_id = "04f00a91"
pdb_path = f"rfd2/out/{config_id}/{config_id}_design_1.pdb"
trb_path = f"rfd2/out/{config_id}/{config_id}_design_1.trb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("design", pdb_path)

with open(trb_path, "rb") as f:
    trb_data = pickle.load(f)
```
