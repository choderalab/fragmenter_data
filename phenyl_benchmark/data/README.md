## Manifest

* `fgroup_R1_wbos.json` - `[wbo, name, SMILES]` for each molecule in fgroup. The R1 bond (or the bond of the phenyl ring
to the R1 functional group) is tagged with a 4 and 5. The list is sorted by WBO in increasing order
* `fgroup_td_job_indices` = `[job_index, [dihedral], name, wbo]` - The QCArchive job index for the selected molecules with 
a bit more information. This is helpful metadata when analyzing QCArchive torsiondrives. 
* `phenyl_set_torsiondrive_inputs.json` - QCArchive torsiondrive inputs
* `qcarchive_torsiondrives.json` - QCArchive torsiondrive results (energies and Lowding-Wiberg bond orders)