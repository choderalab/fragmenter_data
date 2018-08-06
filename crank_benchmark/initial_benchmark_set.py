from fragmenter.workflow_api import workflow
import oenotebook as oenb

# List of molecule to fragment. Both of these molecule only generate one fragment
smiles = ['CCCC', 'c1ccccc1(C(=O)(N))']
# Default crank jobs. Each job gets written out to its own directory with a crank initial state json
crank_jobs = workflow(smiles)
