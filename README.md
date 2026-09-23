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
- [x] Invisible clashes (combining A and B extensions together produces clashing) 
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
python -m rfd2_filtering.filtering --dry-run --bound frac_loop 0.0 0.4 --bound rg_A 0.0 19.0 --bound rg_B 0.0 19.0
# 3.a tip: to list available parameters for filtering, run:
python -m rfd2_filtering.filtering --list-bounds
# 3.b When ready to filter, run without --dry-run flag:
python -m rfd2_filtering.filtering --bound frac_loop 0.0 0.4 --bound rg_A 0 19.0

```

## Step 4 - MPNN
### TLDR:
1. Run MPNN:
```bash
python -m mpnn.write_mpnn_inputs --pdb-dir rfd2/out --passing-csv-path rfd2_filtering/passing.csv
mpnn/run_mpnn.sh /projects/arvid/c_symmetry_design/mpnn/input_files/YYYYmmDD_HHMMSS/
```
2. Run MPNN filter + split
```bash
echo "Yet #TODO"
```


### Rationale
RFD2 output can be sketched like this (pseudo-sequence):
  mb-M-MA--A B--MB-M-ma (1)
  or 
  A--MA-M-mb ma-M-MB--B (2)
If we consider the first case we want to feed the whole thing into MPNN (maybe omitting ligands), but
we do not care for what MPNN predicts about mb on the A chain and ma on the B chain, because in the end,
the two chains will be fused into:
  B--MB-M-MA--A (1)
  or
  A--MA-M-MB--B (2)
We are then left needing no ma nor mb, so they're of no importance.
There are however five, 5, important regions on the rfd2 output left that need to be carefully conveyed
to MPNN with the right instructions for each of them:
- M
- MA
- A
- MB
- B

The real aspects of running this is introduced below: 

### MPNN and indexes
We distinguish the parts of the sequence that matter, and which residues should be redesigned by utilising the --redesigned_residues flag. 
We also use four of the arrays/lists that we get from the trb data filed associated with each output from RFD2: 
- contigmap_hal
- contigmap_ref 
- inpaint_seq
- inpaint_str 

The instructions for each of the five regions are:
1. M
> M is the original input structure, minus any inpaint sequence (MA, MB).
> Notice that we don't care for ma or mb as they can be part of M without interfering with the output, since they will be overwritten by their designed counterparts, MA or MB, respectively.
2. MA + MB, A + B
> If `inpaint_seq[i] == False` the sequence identity of the residue at index *i* is not known.
> It is however possible to discern any inpainted sequence MX from an extension X by consulting inpaint_str:
> If `inpaint_str == True` the residue is MX and otherwise (`inpaint_str == False`) it is an extension X.
> This is nice information to have but it doesn't matter as much for mapping into MPNN as both residue types will get the redesigned_residues label.

The resulting code is run like this:
```bash (ca_rfd)
python -m mpnn.write_mpnn_inputs --pdb-dir rfd2/out --passing-csv-path rfd2_filtering/passing.csv
```


### MPNN filtering
> Q: There's rfd2_filtering, why is filtering of MPNN not it's own submodule?
> A: It's not so complex.



# TODO
- [x] Make a script to easily create (and later access) score distributions.
  - For now `kitten icat *` (in the generated plot folder) or copying the files to local (rsync) is convenient enough.
- [x] Determine if dssp lumps should be part of loop fraction in filtering score.
  *"Yes, Helix, Strand, Loop is conventional, although not a good split for Ordered/Unordered. Decision for now
  is to skip implementing this."*
- [ ] Chain break between N-CO is not caught in the CA-CA distance filter. 
- [ ] Set up to run for all aligned input structures. (`build_rfd2_configs.py`)
- [ ] Set up scoring while running inference. Preferably simultaneously with inference,"on another thread". 
- [ ] Refactor ./input_structures, ./output_structures and ./motif_pipeline into the same subdirectory.

