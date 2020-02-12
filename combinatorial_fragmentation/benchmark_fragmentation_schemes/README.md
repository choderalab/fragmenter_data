## Benchmark different fragmentation schemes on benchmark set

The benchmark set was chosen by finding the molecules and bonds that produce combinatorial fragments that have WBO distributions
with highest distance scores.

To benchmark each scheme:
1. Fragment the selected bonds with scheme
2. For each fragment generate distribution of WBOs from fragments
3. Get distance of distribution to parent
4. Plot against computational cost of fragment (NHeavy)^2.6^

To verify benchmark
1. Besides Omega conformers, generate torsion grid conformers so we get higher energy conformers that are relevant to torsion scans
2. Get distances of 0.03 to parent and Pfizer fragment to parent
3. Compare the distances (This histogram is used in figure 13 in manuscript)
4. By inspection, find functional groups that have strong long distance through-bond effects (shown in figure 14 in manuscript)


## Manifest
### Scripts (in order they were run)
* `calculate_fragment_wbo_dist.py` - script to fragment molecules using different thresholds and generate wbo distributions
* `pfizer_fragmentation.py` - script to generate fragments and their WBO distribution using Pfizer's fragmenation scheme described in https://pubs.acs.org/doi/10.1021/acs.jcim.9b00373
* `torsion_scan_wbos.py` - script to generate conformers from a torsion scan to add some higher energy conformers to WBO distributions to verify the results of
the benchmark.
* `summarize_benchmark.py` - script to generate 2D plots of benchmarking results with fragment computational cost vs distance score
* `compare_scores.py` - script to generate score comparison figure
* `generate_figure.py` - script to generate distributions and fragment visualization that were used to generate manuscript figure.
*`{NAME}/optimal_fragment.py` - script to generate figure 15 in manuscript showing where WBO fragmentation fails.
### scripts that were used for exploration but results were not used for final figures
* `pareto_front.py` - script to generate Pareto front of all fragments from combinatorial fragmentation that have all 1-4 atoms with fragmenter
results plotted in red for visualization (not used to generate figure but helpful to look at).
* `visualize_results.py` - script to visualize results from benchmarking (rudimentary - was not used to generate figures)
### PDF figures
* `jointplot_{}.pdf.format(threshold)` - figures generated with `summarize_benchmark.py`
* `combined-score-differences.pdf` - figure generated with `compare_scores.py`

__note__:
Decisions I made based on intermediate results.
1. While there are several options with fragmenters, I settled on the defaults which are:
    1. Tag functional groups because otherwise fragmenter can generate fragments with weird chemistries
    2. `keep_non_rotor_ring_substituents` set to `False`. There is no need to add these before building out fragment.
    It just leads to larger fragments in general.
    3.  `heuristic` - `path_length`. This led to smaller fragments in general than using the `wbo` heuristic without significantly changing the distance scores.
2. I found a discrepancy with WBO generated on Linux vs Mac for nitrile guanidinium. WBOs on Linux look more correct than Mac so only WBOs that were run on the cluster should be used.

I will need to come back to this to generate SI figures for the different parameters to show that we should use default parameters

How should we handle carbonyl on Pixantrone? It is meta to central bond, but it is on a fused ring so should probably be included.
