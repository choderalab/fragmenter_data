## Combinatorial fragmentation of validation set

Scripts for fragmenting validation set.
All molecules were fragmented along every bond that is not in a ring. All possible
fragments were generated with the following (default) setting:
1. min_rotors: 1
2. max_rotors: rotors in molecule + 1
3. min_heavy_atoms: 4

## Manifest
* `fragment_validation_set.py` - script to generate fragments
* `submit_fragment.py` - script to submit fragment jobs to the queue
* `fragment_validation_set.lsf` - LSF dummy submit script