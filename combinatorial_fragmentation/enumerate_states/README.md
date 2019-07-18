## Enumerate states

Enumerate states for validation set inputs
Since stereo changes do not change WBOs - only one stereoisomer was generated. 
[`OEGetReasonableTautomers`](https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonFunctions/OEGetReasonableTautomers.html#OEProton::OEGetReasonableTautomers)
was used to generate reasonable tautomers. When compared to populations predicted
by chemicalize for pH ~ 7.4, `OEGetReasonableTautomers` is more conservative so we might have 
missed some reasonable states. 

## Manifest
* `enumerate_states.py` - script to enumerate tautomer states at pH ~7.4 for validation set inputs
* `enumerated_states.pdf` - PDF of states
