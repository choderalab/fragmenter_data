## Expand molecule 

Before a molecule gets fragmented, we need to enumerate its protonation
states, tautomers and stereoisomers so that the database covers all torsions
that are likely important. This folder contains the enumeration from 
`fragmenter` with the default settings.

### Notes
* Tautomer enumeration is set to False by default. That is because OpenEye
tautomer enumeration many times gives unreasonable results or enumerates 
resonance structures. 

* Currently all protonation states are enumerated. However, some of them 
are only likely present at very high or low pH. In the future we may want to 
include only states that are likely present in pH ranges that are usually
simulates (6-8)

* For bespoke torsions, it is not always required to expand the molecule if
the user knows which state they want to simulate

### Manifest
* `Enumerate_tautomers.ipynb` - Jupyter notebook to expand molecule
* `kinase_inhibitors.smi` - input smi file for expanding molecule
* `states.smi` - ouput file of expanded states
* `states.csv` - csv file for `csv2xlsx.py` 
* `csv2xlsx.py` - script to visualize results in xlsx file
* `states.xlsx` - visualization of expanded states