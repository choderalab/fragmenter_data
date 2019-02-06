## Torsion scan of ethylene glycol

The global minimum of ethylene glycol is not the trans conformation but the
cis conformation can form intramolecular hydrogen bonding. Here I tested
the wavefront propagation of torsiondrive to find the global minimum.

1. Using only one omega generated starting conformation did not find the global minumum
Associated files:
    * `initial_conformation.xyz` - the initial conformation used for the 1D torsion drive
    * `ethylene_glycol.ipynb` - ipython notebook of submitting the scan to qcfractal 
    * `psi4_final_energies.json` - final energies from QCFractal
    * `psi4_final_molecules.json` - final conformations from QCFractal
    * `workflow.json` - workflows for the scans
    * `visualize_scan.ipynb` - ipython notebook to visualize results
    * `ethylene_glycol_psi4_traj.xyz` - xyz trajectory of final molecules. This was used in pymol to generate images for the movie
    * `movie.mp4` - movie of scan. This movie was generated from pymol images of the conformations

