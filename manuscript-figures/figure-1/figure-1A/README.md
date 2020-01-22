# Generate figure 1 to show how QC torsiondrives scale with molecular size

## Manifest
### scripts
* `td_drives_stats.py` - script to download OFF torsiondrives and save stats and metadata - these sets have pretty small fragmetns,  up to 19 heavy atoms
* `download_ki_opt_data.py` - script to download kinase inhibitor optimization data. This set has larger molecules (23-37 heavy atoms)
but are only optimizations. No torsiondrives

### outputs
* `off_torsions_cpu_time.json` - metadata and cpu time for OFF torsiondrives
* `ki_opt_cput_time.json` - metadata and cpu time for kinase inhibitor optimization dataset
* `opts_per_td.json` - optimization evaluation per torsiondrive and the heavy atoms of the molecules to see if there is a relationship with size
* `gradients_per_opt.json` - gradient evaluations per optimization for OFF torsiondrives and heavy atoms to see if there is a relationship
* `gradients_per_opt_ki.json` - gradient evaluations per optimization for kinase inhibitor optimizations