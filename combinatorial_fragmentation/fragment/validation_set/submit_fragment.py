import os
import subprocess
import pandas as pd

validation_set = pd.read_csv('../../filter/validation_set.csv')

jobs_to_run = validation_set.name

for name in jobs_to_run:
    print(name)
    with open('fragment_validation_set.lsf', 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('JOB_NAME', name)
    bsub_file = os.path.join(os.getcwd(), '{}_fragment.lsf'.format(name))
    with open(bsub_file, 'w') as f:
        f.write(filedata)

    #submit job
    stdin_file = open(bsub_file, 'r')
    subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
    stdin_file.close()

