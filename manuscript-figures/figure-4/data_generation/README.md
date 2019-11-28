## Generate distributions of WBO for a set of kinase inhibitors

## Manifest
* `archive/` - initial hf3c calculations. The calculations were rerun on QCArchive
* `data/` - mol2 and oeb files of conformation used to generate AM1 inputs and QCArchive inputs. oeb files have AM1 WBOs
* `kinase_inhibitors.smi` - kinase inhibitors used for this study. This file does not include inhibitors with large ring sytems
such as Everolimus, Lorlatinib, Midostaurin, Sirolimus, Temsirolimus, Alectinib
* `visualize_kinase_inhibitors.py` - script to visualize kinase inhibitors
* `kinase_inhibitors.pdf` - visualization of kinase inhibitors used to generate data
* `01_generate.py` - script to generate qcarchive geometry optimizaiton input. to run `python 01_generate.py -i kinase_inhibotors.smi`
* `02_create_optimization_input.py` - script to generate QCArchive OptimizationDataset

__note__: The scripts to generate and submit the conformers for calculations on QCArchive are [here](https://github.com/openforcefield/qca-dataset-submission/pull/69)
* `calculate_am1_wbo.py` - script to calculate AM1 WBOs
* `submit.py` - script to submit AM1 calculation on cluster