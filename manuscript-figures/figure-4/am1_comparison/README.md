## Compare AM1 WBO distributions and correlation to HF3C

The original WBO distributions were calculated using HF3C to make sense of AM1 ELF10 WBOs. When looking at these
distributions I found that the variance of WBO was higher for conjugated bonds and bonds in conjugated systems have WBOs
that are correlated to each other. In this folder I wanted to see if this trend is reproducible with AM1 WBOs for each
conformation. I tested Imatinib and Gefitinib and found that the trend is somewhat reproducible but does look different.
Conjugated bonds have higher variance but the variance in general is higher than in HF3C. The correlations of WBOs in
conjugated systems are not as strong as HF3C

## Manifest
* `calculate_am1_wbo.py` - script to calculate wbos
* `Gefitinib_for_am1_comparison.oeb` - gefitinib oemols with WBOs (output from running `python calculate_am1_wbo.py -n Gefitinib`)
* `Imatinib_for_am1_comparison.oeb` - Imatinib oemols with WBOs (output from running `python calculate_am1_wbo.py -n Imatinib`)
