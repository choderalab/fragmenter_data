import os
import subprocess
import glob
import itertools

jobs_to_run = glob.glob('selected/*/')
parameters = {'threshold': [0.5, 0.1, 0.05, 0.01, 0.005, 0.001],
              'path': ['path_length', 'wbo'],
              'functional_groups': ['None', 'False'],
              'keep_non_rotor': ['True', 'False']}
for dir in jobs_to_run:
    name = dir.split('/')[-2]
    print(name)

    with open('fragment_scheme_score.lsf', 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('JOB_NAME', name)
    i = 0
    for t in parameters['threshold']:
        for p in parameters['path']:
            for fg in parameters['functional_groups']:
                for nr in parameters['keep_non_rotor']:
                    with open('fragment_scheme_score.lsf', 'r') as f:
                        filedata = f.read()
                    filedata = filedata.replace('JOB_NAME', name)
                    filedata = filedata.replace('THRESHOLD', str(t))
                    filedata = filedata.replace('PATH', p)
                    filedata = filedata.replace('FUNCTIONAL_GROUPS', fg)
                    filedata = filedata.replace('NONROTORS', nr)

                    bsub_file = os.path.join(os.getcwd(), '{}_fragment_scheme_score_{}.lsf'.format(name, str(i)))
                    with open(bsub_file, 'w') as f:
                        f.write(filedata)
                    i +=1

                    #submit job
                    stdin_file = open(bsub_file, 'r')
                    subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
                    stdin_file.close()
