import fragmenter
import pandas as pd
from qcfractal import interface as ptl
import json

validation_set = pd.read_csv('../filter/validation_set.csv')
client = ptl.FractalClient('https://localhost:7777/', verify=False)

workflow = {
    "validation_set": {
        "fragmenter":{
          "enumerate_states": {
            "version": fragmenter.__version__,
            "options": {
              "protonation": True,
              "tautomers": False,
              "stereoisomers": False,
              "max_states": 10,
              "level": 0,
              "reasonable": True,
              "carbon_hybridization": True,
              "suppress_hydrogen": True,
              "mapped": True
            }
          },
          "enumerate_fragments": {
            "version": fragmenter.__version__,
            "scheme": "combinatorial",
            "functional_groups": False,
            "options":{},
          },
        "torsiondrive_input": {},
        "torsiondrive_static_options": {
        "optimization_spec": {
          "program": "geometric",
          "keywords": {
            "coordsys": "tric"
          }
        },
        "qc_spec": {
          "driver": "gradient",
          "method": "UFF",
          "basis": "",
          "keywords": None,
          "program": "rdkit"
        },
        "keywords": {}
      },
      "optimization_static_options":{
          "program": "geometric",
          "keywords": {
            "coordsys": "tric"
          },
        "qc_spec":{
           "driver": "gradient",
          "method": "UFF",
          "basis": "",
          "keywords": None,
          "program": "rdkit"
        }
      }
    }
  }
}


with open('../validation_set_workflow.json', 'w') as f:
    json.dump(workflow, f, sort_keys=True, indent=2)

fragmenter_workflow = fragmenter.workflow_api.WorkFlow(workflow_id='validation_set', client=client,
                                                       workflow_json='../validation_set_workflow.json')

for smiles, name in zip(validation_set.smiles, validation_set.name):
    fragmenter_workflow.enumerate_states(smiles, title=name, json_filename='states_{}.json'.format(name))