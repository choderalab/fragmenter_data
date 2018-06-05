## Mayer and Wiberg fractional bond order

Configurations were generated with Omega. 
Wiberg Lowden and Mayer indices were calculated at hf3c/def2-svp level of theory
for all configurations. 

### Manifest
* `psi4_input.py` - script that generated psi4 inputs and submitted jobs on cluster
* `bond_order_bsub.lsf` - dummy lsf script to submit jobs
* `run_psi4.py` - script to launch psi4 and save output
* `psi4_bond_order.ipynb` - notebook used to generate plots. 

* `*_bond_orders.pdf` - pdf file of bond order results for FDA approved kinase inhibitors. 
 The title of each histogram corresponds bond atom indices which corresponds to the numbers
 on the molecule image. The histogram is over the configurations of the molecule
* `*_bond_orders_no_ring.pdf` - This file only plots the bonds of interest. Since we are
not fragmenting rings their bond order is not that important here. The black line corresponds
to OpenEye Wiberg bond order. The red line corresponds to 1.2 which is our current threshold.
The interesting bonds are the ones where the threshold is between OpenEye and/or psi4 bond orders.
The bonds that seem to give results like that are bonds that include nitrogen and halogens. 