## Benchmark different fragmentation schemes on benchmark set

The benchmark set was chosen by finding the molecules and bonds that produce combinatorial fragments that have WBO distributions
with highest distance scores.

To benchmark each scheme:
1. Fragment the selected bonds with scheme
2. For each fragment generate distribution of WBOs from fragments
3. Get distance of distribution to parent

Figures to do:
1. 2D plot of size vs distance for all parameters in 0.1, 0.05, 0.01, 0.005, 0.001
2. Find the one 2D that looks the best and generate plots separating out different parameters (functional groups, rotors, huerisitc (6 2D plots))

For paper
Figure 1 - all molecules in benchmark set (supplementary)
Figure 2 - example of combinatorial fragmentation of a molecule that has good binning. Pareto front
Figure 3 - summary of benchmark for current scheme with different parameters. Supplementary separate out different parameters.
Figure 4 - show an example of figure 3 in benchmark (nice example:  (Abemaciclib (34, 16), Bosutinib (30, 10), Cefepime_1 (19, 25), Ceftazidime_0 (11, 7) (25, 35)?
Cibmimetnib_0 (13, 22) Dacomitinib_0 (28, 17), Aceclofenac_0 (14, 21),
Also - add Bemis Murco fragmentation, Pfizer fragmentation

Bonds in molecules where you inevitably get the parent back even though there is a more optimal fragment, this scheme is not finding it
0.1 - Fostamatinib_0 (29, 14), Pemetrexed_0 (13, 25)
0.05 - Abemaciclib_0 (14, 34), Bosutinib_0 (11, 30)

Bonds in molecule where score is too high for 0.05.
Hydroxychloroquin_0 (8, 20) - not finding optimal.
Netarsudil_0 (30, 23)
Proguanil_4 (7, 13)

Where Pfizer fails
1. Menadiol_diphosphate - need phosphate group on other side.
2. Gemifloxacin (8, 22) no one gets. Pfizer does worse in the score but fragmenter does worse in the size.
3. Phenformin_3 (13, 7) - Pfizer drops a positively charged nitrogen
4. Proguanil_0 (14, 7) Pfizer drops a positively charged nitrogen (ridge plot is better for figure)
5. Tiludronic_acid_0 (16, 7) - loss of chlorine seems to matter here.


## Manifest
* `calculate_fragment_wbo_dist.py` - script to fragment molecules using different thresholds and generate wbo distributions
* `pfizer_fragmentation.py` - script to generate fragments and their WBO distribution using Pfizer's fragmenation scheme described in https://pubs.acs.org/doi/10.1021/acs.jcim.9b00373
* `summarize_benchmark.py` - script to generate 2D plots of benchmarking results with fragment computational cost vs distance score
* `pareto_front.py` - script to generate Pareto front of all fragments from combinatorial fragmentation that have all 1-4 atoms with fragmenter
results plotted in red for visualization.
* `visualize_results.py` - script to visualize results from benchmarking
* `torsion_scan_wbos.py` - script to generate conformers from a torsion scan to add some higher energy conformers to WBO distributions to verify the results of
the benchmark.
* `jointplot_{}.pdf.format(threshold)` - figures generated with `summarize_benchmark.py`

__note__:

Files with `fixed` in them denote that they were generated with updates to `fragmenter` because some functional groups were
not in the yaml file so they were fragmented. Also, for the second time around, the parameters for fragmentation were set to the default
values because the previous run showed that those had the best results in general. The default parameters are:
1. `functional_groups`: None - use `fragmenter`s internal yaml file for functional groups not to fragment. If you don't tag these groups you can end up with weird fragments
2. `keep_non_rotor_ring_substituents` - `False`. There is no need to add these before building out fragment. It just leads to larger fragments in general
3. `huerisitc` - `path_length`. This led to smaller fragments than using the `wbo` heuristic.
