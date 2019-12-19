## Benchmark different fragmentation schemes on benchmark set

The benchmark set was chosen by finding the molecules and bonds that produce combinatorial fragments that have WBO distributions
with highest distance scores.

To benchmark each scheme:
1. Fragment the selected bonds with scheme
2. For each fragment generate distribution of WBOs from fragments
3. Get distance of distribution to parent

## Manifest
* `calculate_fragment_wbo_dist.py` - script to generate wbo distributions
