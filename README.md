# c_symmetry_design

## Step 0 - Align monomers
Here we use the **Center Of Mass** (COM) of the protein as a whole and the COM of the *Ligand* to standardise the monomer so that we can push it through the pipeline. This essentially carried out by placing the protein's COM at origin and pointing the ligand (more mathematically a vector from COM_Protein to COM_Ligand is created and aligned) to the positive x-direction. Because one might want to tweak this (e.g. assessing placement in PyMOL), and doing this over and over is a computational cost, the structures are saved in `c_symmetry_design/aligned/`. 

## Step 1 - Generating ppi motifs
Essentially: `python -m motif_pipeline.generate_pdbs`

## Step 2 - Running inference (**RFD2**)
1. Generate configs for RFD2. `python -m rfd2.build_rfd2_configs`
2. Run RFD2. `rfd2/run_all.sh`

## Step 3 - Filtering output
Current scores:
- [x] Clash
- [x] Backbone continuity
- [x] Secondary structure
- [x] Normalised Radius of Gyration
- [x] Center of Mass Shift
- [ ] Interface 
- - [ ] Intertwining
- - [ ] 

1. Check backbone integrity (backbone and clash)
`python -m rfd2_filtering.run_scoring --out-dir rfd2/out --scores-path rfd2_filtering/scores.csv`
2. Secondary Structure using DSSP. 
`python -m rfd2_filtering.run_dssp --scores-path rfd2_filtering/scores.csv --dssp-scores-path rfd2_filtering/dssp_scores.csv`


# TODO
- [ ] Give some sort of life signal when running rfd2_filtering.run_scoring
- [ ] Make a script to easily create (and later access) score distributions.
- [ ] Check clashing of the combined extensions (Align and check for clashes outside original contig)
- [ ] 
