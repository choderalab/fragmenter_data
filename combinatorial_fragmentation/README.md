## Exhuastive fragmentation

### Generate all possible fragments for a set of molecules
We generated this set to investigate how a bond's chemical environment changes with different remote chemical changes.

Steps taken to generate this set:
1. Filter DrugBank. Scripts used are in the `filter/` folder
2. Generate reasonable tautomers for filtered molecules. Scripts used for this step are in the `enumerate_states/` folder. Resulting JSON files are on Lilac.
3. Exhaustive fragmentation. Scripts for this step are in the `fragment/` folder. Resulting fragments are in JSON files on Lilac.
4. Generate conformers and calculate their WBOs. scripts are in the `fragment_bond_orders/` folder. All resulting JSON
files are on Lilac.

### Choose benchmark set.
The goal here is to find molecules that are challenging to fragment by finding fragments where remote
chemical changes cause large changes to the WBO distributions of the central bonds. To do that, each fragment is scored
by its bonds' WBO distribution distance to parent WBO (using the `mmd_x_xsqred` function in the `score_fragments.py` script in `rank_fragments`).
Distance scores were sorted and the bonds with the top 100 scores were selected for the benchmark score.
__note__: Not all bonds in a selected molecules is part of the benchmark set, only the bonds that had distance scores over 0.001 (should be 0.03 with sqrt)
The initial selection had a bug where the sqrt was not taken for the mmd score so the distances are in general lower than they should be.
It was fixed later for the final benchmark of different fragmentation scheme. The selection might be slightly different if this was fixed
before the final 100 moelecules were chosen, but the difference does not change the conclusions of the study because the distances used
for later analysis and the actual benchmark did have the sqrt taken.

* `rank_fragments/` - folder containing scripts to select 100 molecules for the benchmark set

### Benchmark fragmentation schemes
Generate fragments for the benchmark set using different fragmenation schemes.
Schemes that were tested:
1. Scheme described by the Pfizer group
2. WBO fragmentation with different thresholds and parameters.

* `benchmark_fragmentation_schemes/` - scripts to generate fragments using different schemes and analyze results.

__note__:
A good experiment to do with this set is run QC torsion scans of 2 fragments for each molecule
1. Smallest fragment with greatest distance to parent
2. Smallest fragment with largest distance to parent.
This will provide the best evidence for the way we are evaluating fragmentation schemes. I just didn't have time to do it -
reviewers might ask for it though and OFF should do it regardless.

This directory also includes some scripts and data for the original experiment on kinase ihibitors and mini drugbank that
then inspired me to look at the entire DrugBank to find more interesting chemistries that will have long range effects.
