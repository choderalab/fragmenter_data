import os
import subprocess
import glob
import itertools

jobs_to_run = glob.glob('../enumerate_states/validation_set/*_states.json')

for dir in jobs_to_run:
    name = dir.split('/')[-1]
    name = name[:-12]
    print(name)
    with open('fragment.lsf', 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('JOB_NAME', name)

    bsub_file = os.path.join(os.getcwd(), '{}_fragment.lsf'.format(name))
    with open(bsub_file, 'w') as f:
        f.write(filedata)

    stdin_file = open(bsub_file, 'r')
    subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
    stdin_file.close()
