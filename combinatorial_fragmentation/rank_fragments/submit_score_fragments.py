import os
import subprocess
import glob

jobs_to_run = glob.glob('../fragment_bond_orders/validation_set/*/')

for dir in jobs_to_run:
    name = dir.split('/')[-2]
    print(name)

    with open('score_fragments.lsf', 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('JOB_NAME', name)
    bsub_file = os.path.join(os.getcwd(), '{}_score_fragments.lsf'.format(name))
    with open(bsub_file, 'w') as f:
        f.write(filedata)

    #submit job
    stdin_file = open(bsub_file, 'r')
    subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
    stdin_file.close()
