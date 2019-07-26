# Data for the development of `fragmenter`

This repository contains scripts, data and figures to explore the science behind `fragmenter`


## Manifest
* `bond_order_variance/` - scripts and figures exploring Wiberg bond order variance, how it is conformational dependant 
and what that tells us about chemical environments
* `combinatorial_fragmentation/` - scripts and figures for generating the benchmark validation set for `fragmenter`
* `conjugation/`- scripts and figures exploring Wiberg bond order correlations w.r.t. conformation and how 
that can inform us on conjugated systems. Comparisons of different conjugation models.
* `phenyl_benchmark/` - scripts and figures for the phenyl set and biphenyl motivational example. This set explores
how Wiber bond orders change with changes in chemical environment of remote substituents at least 2 bonds away. The set
was initially inspired by the [Hammet equation](https://en.wikipedia.org/wiki/Hammett_equation)
* `frag_opt/` - optimization of `fragmenter` by looking at different paths it can take when growing out a fragment