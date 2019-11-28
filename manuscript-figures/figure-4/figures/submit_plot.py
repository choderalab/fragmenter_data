import fragmenter
import subprocess

oemols = fragmenter.chemi.file_to_oemols('../data_generation/kinase_inhibitors.smi')
for mol in oemols:
    name = mol.GetTitle()
    print(name)
    try:
        subprocess.call("python generate_am1_figures.py -n {}".format(name), shell=True)
    except:
        print('plotting failed for {}'.format(name))
