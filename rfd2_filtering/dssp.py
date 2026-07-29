import os
import tempfile
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.DSSP import DSSP

def fix_tmp_file(file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("my_protein", file)

    # 1. Create and write the file, then close it completely
    with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w+', delete=False) as tmp:
        tmp_path = tmp.name
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp)  # Writes to file handle

    # At this point, exiting the 'with' block flushes and closes 'tmp' on disk

    try:
        # 2. Parse and execute DSSP on the closed file path
        fixed_structure = parser.get_structure("my_fixed_protein", tmp_path)
        model = fixed_structure[0]

        dssp = DSSP(model, tmp_path, dssp='mkdssp')

        for key in dssp.keys():
            residue_id = key
            sec_struct = dssp[key][2] 
            print(f"Residue {residue_id}: SS= {sec_struct}")
 
    finally:
        # 3. Clean up disk entry
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

if __name__ == "__main__":
    local_test_pdb = "/projects/arvid/c_symmetry_design/rfd2/out/0001adcf/design_5.pdb"
    fix_tmp_file(local_test_pdb)
