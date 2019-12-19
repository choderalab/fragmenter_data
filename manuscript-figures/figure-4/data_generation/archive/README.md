## Generating data to analyze conformational dependency of WBOs

__Note__:
The data was generated before QCArchive was available.

The initial reason for running these calculations was to test how close OpenEye ELF10 conformation independent Wiberg
Bond orders are to bond orders calculated with higher level of theory. Conformations were generated with Omega, optimized
with hf3c and WBOs were calculated using Löwdin normalization. The results ended up showing that Wiberg bond orders
calculated with hf3c has higher variance for conjugated bonds and that bond orders in conjugated systems are correlated
with each other. When the same process was later tried with OpenEye AM1 WBO, I was able to reproduce the result of higher
variance for conjugated bond but the conjugated bonds did not show such strong correlations
as the hf3c Wiberg-Löwdin (see AM1_comparison)

## Manifest
* `clinical_kinase_inhibitors_tagged.smi` - tagges SMILES of kinase inhibitors used in this study. Note: These tags
were generated before `cmiles` used canonical atom order so they are ordered differently than current `cmiles` mapped SMILES
* `run_psi4.py` - script that reads psi4 input files and runs psi4 with options specified in psi4 input files
* `psi4_input.py` - script to generate psi4 input files. Omega was used to generate conformers
* `geom_opt_bsub.lsf` - dummy lsf bsub file to run geometry optimization jobs
* `psi4_bo_input.py` - reads geometry output files for optimized geometry, generates psi4 input file with optimized geometry and
options to calculate Wiberg-Löwding and Mayer bond orders
* `bo_bsub.lsf` - dummy lsf bsub script to run bond order calculations

