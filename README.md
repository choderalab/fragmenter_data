# Data for the development of `fragmenter`

This repository contains scripts, data and figures to explore the science behind `fragmenter`

## Manifest
* `fragmenter-manuscript-figures/` - scripts and data used for the final figures that made it into the first manuscript (https://chayast.github.io/frag-manuscript/). 
* `wbo-manuscript-figures/` - scripts and data used for the final figures that made it into the second manuscript (https://chayast.github.io/wbo-manuscript/) (WIP)
* `combinatorial_fragmentation/` - scripts and figures for generating the benchmark exhaustive fragmentation validation set for `fragmenter` (used for first manuscript). 
* `phenyl_benchmark/` - scripts and figures for the phenyl set and biphenyl motivational example. This set explores
how Wiber bond orders change with changes in chemical environment of remote substituents at least 2 bonds away. The set
was initially inspired by the [Hammet equation](https://en.wikipedia.org/wiki/Hammett_equation) (used for second manuscript). 

### Other folders containing scripts and analysis that did not make it into the paper.
Most of these inspired better experiments that did end up in the paper.
* `bond_order_variance/` - scripts and figures exploring Wiberg bond order variance, how it is conformational dependant 
and what that tells us about chemical environments
* `conjugation/`- scripts and figures exploring Wiberg bond order correlations w.r.t. conformation and how 
that can inform us on conjugated systems. Comparisons of different conjugation models.
* `frag_opt/` - optimization of `fragmenter` by looking at different paths it can take when growing out a fragment
