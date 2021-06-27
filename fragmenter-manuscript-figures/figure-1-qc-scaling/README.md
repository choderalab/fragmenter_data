# Generate figure 1 to show how QC torsiondrives scale with molecular size

## Manifest
### scripts
* `td_drives_stats.py` - script to download OFF torsiondrives and save stats and metadata - these sets have pretty small fragmetns,  up to 19 heavy atoms

__note__:
It took forever to download the data so I only managed to download up to molecule 689 in the OpenFF Group 1 Torsions.

* `download_ki_opt_data.py` - script to download kinase inhibitor optimization data. This set has larger molecules (23-37 heavy atoms)
but are only optimizations. No torsiondrives
* `dft_scaling.py` - script to generate figure showing DFT scaling, distributions of opts/torsiondrive and gradient evaluations / optimization

### outputs
* `off_torsions_cpu_time.json` - metadata and cpu time for OFF torsiondrives
* `ki_opt_cput_time.json` - metadata and cpu time for kinase inhibitor optimization dataset
* `opts_per_td.json` - optimization evaluation per torsiondrive and the heavy atoms of the molecules to see if there is a relationship with size
* `gradients_per_opt.json` - gradient evaluations per optimization for OFF torsiondrives and heavy atoms to see if there is a relationship
* `gradients_per_opt_ki.json` - gradient evaluations per optimization for kinase inhibitor optimizations

### Figures
* `B3LYP_scaling.pdf` - figure showing DFT scaling for one CPU and distribution of druglike molecule size by heavy atoms
* `SI_Intel_scaling.pdf` - SI figure showing DFT scaling for other Intel Xeon processors
* `SI_AMD_scaling.pdf` - SI figure showing DFT scaling for other AMD processors
* `opts_per_td.pdf` - boxplot of optimizaitons per torsion drive by molecular size
* `gradients_per_opts.pdf` - boxplot of gradient evaluations per optimization by molecular size