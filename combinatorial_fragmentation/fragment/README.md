## Combinatorial Fragmentation

Combinatorial fragmentation of molecules. Molecules are fragmented along
every bond besides bonds in rings and selected functional groups. Then 
all possible connected fragments are generated combinatorial between 
specified minimum and maximum rotatable bonds

## Manifest
* `kinase_inhibitors/` combinatorial fragments from the filtered kinase inhbitor set  
For this set, the following functional groups were not fragmented:
1. amide
2. Triflouro methyl
3. Sulfonamide
4. Sulfone
5. Ester
6. Phosphate
7. P(=O)(C)(C) 
* `validation_set/` combinatorial fragmentation from the filtered drugbank set  
For this set, all functional groups were fragmented besides rings. 
* `combinatorial_fragmentation.ipynb` - notebook that was used to fragment kinase inhibitor set