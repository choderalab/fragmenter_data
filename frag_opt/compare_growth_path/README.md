## Compare different paths of growing out fragment

More than one path exists to grow out fragment if WBO changes more than the threshold. 
Here I looked at the stability of two heuristics:

1. Shortest path length 
Grow out the bond that is closest to the central bond the fragment is being built for. 
If more than one point to grow exists, choose the bond that has the higher Wiberg Bond Order
The rationale here is that closer atoms will effect the central bond more.
2. Highest Wiberg Bond Order
Grow out the fragment with the bond that has the highest Wiberg Bond Order. The rationale is 
that a bond with greater WBO is probably more conjugated so will contribute or withdraw more 
electronic density from the central bond

### Results
Before running these experiments I thought that the WBO of the central bond will
stabilize at some fragment size and one path will converge faster than the other. Then we
will run QM torsion scans on selected fragments and see if the QC scans stabilize. Turns
out the WBO can oscillate and larger fragments are not necessarily better than smaller ones. 

The filtered kinase inhibitor set was used. The smallest fragment contains all 1-4 atoms around
the central bond. Then one bond is added at a time depending on the the path chosen. 

For most bonds - the WBO did not change that much. To generate the set of fragments to run 
torsion scans on the fragments were filtered to the central bonds that exhibited changes in WBO > 0.1

However, since most of these bonds were between rings (either one or two of the atoms in the 
rotatable bond is in a ring), the scans had difficulty converging. Therefore we are 
currently on running Idelalisib because it is the least restricted rotatable bond. 

### Manifest
For each kinase inhibitor:
* `kinase_inhibitor/kinase_inhibitor_fragment_growth.pdf` - growth curves
* `kinase_inhibitor/kinase_inhibitor_bond_a1_a2.pdf` - fragments corresponding to data points on growth curve
* `kinase_inhibitor/kinase_inhibitor_to_drive.json` - mapped SMILES and some provenance for fragments in growth curve
* `growth_plot.py` - script to generate growth curves
* `generate_stability_benchmark.ipynb` - notebook that filters the fragments around the bonds that exhibit
changes in WBO . 0.1
* `Idelelalisib_frags.json` - Idelelalisib fragments to run
* `kinase_inhbitor_bond_a1_a2_to_run.pdf` - the other fragments that should be run 