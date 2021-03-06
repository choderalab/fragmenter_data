## Rank fragments

This directory contains scripts and figures for selecting the final validation / benchmark set
and working out the score of different fragmentation schemes.  
To select the set, each fragment was scored with the maximum mean discrepancy of its WBO
distribution to the parent WBO distribution using the mean and squared mean distance.
Then, the molecules were sorted by their mmd score and the top 100 scoring fragments
were selected. 

## Manifest
* `score_fragments.py` - script to score fragments and generate ridge plots of the results colored 
with the score. This script only scores fragments that have all 1-5 atoms around central bond
* `rescore_fragments.py` - script to score fragments. This script scores fragments that has bond of interest, regardless of how many
atoms is has around the bond. These results were NOT used to select final benchmarking sets. It was just used on the selected set just to see if the scores
should be normalized. Pretty worthless and should probably go in archive.
* `select_valdiation_set.py` - script to sort the scores and select to 100 molecules
* `selected_validation_set.pdf` - visualization of the selected molecules with the bond with the 
highest mmd score highlighted
* `pareto_front.py` - script to generate the Pareto front for the top 100 molecules