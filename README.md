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
- [ ] Invisible clashes (combining A and B extensions together produces clashing) 
- [ ] Interface 
- - [ ] Intertwining
- - [ ] 

1. Run scoring (skipping dssp for clashing/broken backbones)
Here's my current workflow to check output from RFD2:

```bash
# 1. Run Scoring to measure properties of generated backbones
python -m rfd2_filtering.run_scoring --out-dir rfd2/out --scores-path rfd2_filtering/scores.csv
# 2. Check the outputs by plotting and viewing (make sure to `kitten ssh lama` to use icat in kitty)
python -m rfd2_filtering.plot_scores
kitten icat rfd2_filtering/plots/YYYYMMDD.../*
# 3. Run filtering with dry-run to check passing backbones
python -m rfd2_filtering.filtering --dry-run --bound frac_loop 0.0 0.4 --bound rg_A 0 19.0
# 3.a tip: to list available parameters for filtering, run:
python -m rfd2_filtering.filtering --list-bounds
# 3.b When ready to filter, run without --dry-run flag:
python -m rfd2_filtering.filtering --bound frac_loop 0.0 0.4 --bound rg_A 0 19.0


```



# TODO
- [x] Make a script to easily create (and later access) score distributions.
- - For now `kitten icat *` (in the generated plot folder) or copying the files to local (rsync) is convenient enough.
- [x] Safe delete standalone dssp script in rfd2_filtering.
- [ ] Determine if dssp classified loops should be part of loop fraction in filtering score.
- [ ] Check clashing of the combined extensions (Align and check for clashes outside original contig)
- [ ] Chain break between N-CO is not caught in the CA-CA distance filter. 
- [ ] Invisible clashes (combining A and B extensions together produces clashing) 
- [ ] Set up to run for all aligned input structures. (`build_rfd2_configs.py`)
- [ ] Set up scoring while running inference. Preferably simultaneously with inference,"on another thread". 
