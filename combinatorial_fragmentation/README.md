## Brute force fragmentation 
### All possible fragments from parent molecule

To find how well WBO correlates with change in torsion profiles, set of 
kinase inhibitors are fragmented combinatorially to find all possible fragments. 

Full combinatorial fragmentation and running torsion scans on all is too expensive.
For the kinase inhibitor set it generated 24,000 fragments. To decrease cost:  

1. Filter kinase inhibitor 
    1. Less than 40 atoms 
    2. At least 4 rotatable bond and at most 8 rotatable bonds
    3. Ring systems up to 14 atoms
    4. Remove carbon chains with more than 4 carbons
    (For the first try I removed Trametinib from the list because of the amount of fragments it produced with combinatorial fragmentation)
    
1. Create list of functional groups not to fragment (lives in yam file)
    1. amide
    2. Triflouro methyl
    3. Sulfonamide
    4. Sulfone
    5. Ester
    6. Phosphate
    7. P(=O)(C)(C) 
2. Combinatorially fragmment without fragmenting functional group in list
3. Calculate WBO for all fragments
4. For each rotatable bond in parent molecule:
    1. Find pair of bonds that are closest in molecular weight but furthest in WBO
    2. Run 1D torsion profile of those bonds