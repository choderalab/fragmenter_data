import json
import os
import subprocess
from fragmenter import utils
from openeye import oechem
from openmoltools import openeye
import glob
from cclib.parser import Psi

output_path = os.path.abspath('.')
input_path = os.path.abspath('../clinical_kinase_inhibitors_tagged.smi')
# Load mols from .smi
mollist = utils.to_oemol(input_path, title='title', verbose=False)

json_data = {}
json_data["driver"] = "property"
json_data["kwargs"] = {"properties": ["WIBERG_LOWDIN_INDICES", "MAYER_INDICES"]}
json_data["method"] = "hf3c"
#json_data["options"] = {"BASIS": "def2-svp"}
json_data["return_output"] = "True"

for mol in mollist:
    name = mol.GetTitle()
    print('\nGenerating input for {}'.format(name))
    directory_path = os.path.join(output_path, name)
    try:
        os.mkdir(directory_path)
    except FileExistsError:
        print("Overwriting existing directory")
        pass
    # Load all psi4 output files
    psi_out_txt = glob.glob('{}/*output.txt'.format(directory_path))
    json_out =open("{}/{}_0.output.json".format(directory_path, name), 'r')
    data  = json.load(json_out)
    json_out.close()
    json_data["tagged_smiles"] = data["tagged_smiles"]
   
    for txt in psi_out_txt:
        bo_infile = txt.replace('output.txt', 'input.bo.json')
        log = Psi(txt)
        psi_data = log.parse()
        # Check if geom converged
        try:
            psi_data.optdone
        except AttributeError:
            print('Geometry did not converge for {}'.format(txt))
            continue
        
        json_data['molecule'] = psi_data.writexyz()[15:]
        json_path = os.path.join(name, bo_infile)
        print(json_path)
        with open(json_path, 'w') as outfile:
            json.dump(json_data, outfile, indent=4, sort_keys=True)
        bo_infile = bo_infile.split('/')[-1]
        print(bo_infile)
        # replace command on psub script
        with open('bo_bsub.lsf', 'r') as file:
            filedata = file.read()
        i = txt.split('/')[-1].split('.')[0].split('_')[-1]
        JOB_NAME = '{}_{}'.format(name, str(i))
        filedata = filedata.replace('JOB_NAME', '{}'.format(JOB_NAME))
        
        filedata = filedata.replace('INPUTFILE', '{}'.format(bo_infile))
        outfile = '{}_{}.output'.format(name, str(i))
        filedata = filedata.replace('OUTPUTFILE', '{}'.format(outfile))

        bsub_file = os.path.join(os.getcwd(), '{}/{}_{}_bo.lfs'.format(name, name, str(i)))
        with open(bsub_file, 'w') as file:
            file.write(filedata)
        #change directory and submit job
        os.chdir(os.path.join(os.getcwd(), name))
        print(os.getcwd())
        stdin_file = open('{}_{}_bo.lfs'.format(name, str(i)), 'r')
        subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
        stdin_file.close()
        os.chdir("..")

# import numpy as np
# times = np.asarray(atom_map_time)
# np.save('ss_time', times)
# import matplotlib.pyplot as plt
# plt.hist(atom_map_time)
# plt.xlabel('Time (seconds)')
# plt.title("distribution of substructure search time")
# plt.savefig('ss_time.png')





