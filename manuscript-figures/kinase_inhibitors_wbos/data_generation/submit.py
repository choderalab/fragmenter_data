import fragmenter
import subprocess
import os

oemols = fragmenter.chemi.file_to_oemols('kinase_inhibitors.smi')
for mol in oemols:
    name = mol.GetTitle()
    print(name)
    # replace command on psub script
    with open('calc_am1_bsub.lsf', 'r') as file:
        filedata = file.read()
        filedata = filedata.replace('JOB_NAME', '{}'.format(name))
        bsub_file = os.path.join(os.getcwd(), '{}_calc_am1_bsub.lfs'.format(name))
        with open(bsub_file, 'w') as file:
            file.write(filedata)

        stdin_file = open('{}_calc_am1_bsub.lfs'.format(name), 'r')
        subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
        stdin_file.close()

