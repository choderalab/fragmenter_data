## Figure 2

This folder contains all scripts, data (or how the data was accessed) for figure 2.

__note__:
Used DrugBank Version 5.1.3 (release date: 2019-04-02) Downloaded on 2019-06-07. Extraction and properties was
done in `figure-1/extract_small_mols.py`

## Manifest
### Scripts
* `calculate_elf10_am1_wbos.py` - script to filter drug bank molecules for molecules with tiple bonds and generate wbos those molecules and kinase inhibitors
*  `generate_distributions.py` - script to generate distributions of ELF10 WBOs

### Input files
* `kinase_inhibitors.smi` - kinase inhibitors SMILES
* `drugbank_small_mols.csv` - DrugBank small molecules (Downloaded and filtered on 2019-06-07)

### Output files
* `druglike_wbos.oeb` - drug like oe molecules with wbos
* `druglike_wbos.smi` - SMILES for molecules used to generate WBO distributions
* `druglike_mols_for_wbo_dists.pdf` - PDF of molecules use (use in SI)
