# c_symmetry_design

## Step 0 - Align monomers
Here we use the **Center Of Mass** (COM) of the protein as a whole and the COM of the *Ligand* to standardise the monomer so that we can push it through the pipeline. This essentially carried out by placing the protein's COM at origin and pointing the ligand (more mathematically a vector from COM_Protein to COM_Ligand is created and aligned) to the positive x-direction. Because one might want to tweak this (e.g. assessing placement in PyMOL), and doing this over and over is a computational cost, the structures are saved in `c_symmetry_design/aligned/`. 

## Step 1 - Generating ppi motifs
Essentially: `python -m motif_pipeline.generate_pdbs`

## Step 2 - Running inference (**RFD2**)
1. Generate configs for RFD2. `python -m rfd2.build_rfd2_configs`
2. Run RFD2. `rfd2/run_all.sh`

# TODO
- Score interface by secondary structures - DSSP
- Filter interface by clashes


