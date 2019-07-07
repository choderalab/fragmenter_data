import os
import subprocess
import glob
# jobs_to_run = ['Larotrectinib', 'Lenvatinib', 'Nilotinib', 'Palbociclib', 'Pazopanib', 'Ponatinib', 'Ribociclib',
#                'Sunitinib', 'Tofacitinib', 'Vemurafenib']

jobs_to_run = glob.glob('DrugBank*')

name_seen = []
for file in jobs_to_run:
    name = file.split('_')[1].split('_')[0]
    name = 'DrugBank' + '_' +  name
    if name not in name_seen:
        print(name)
        name_seen.append(name)
        with open('oe_wbo_parent.lsf', 'r') as f:
            filedata = f.read()
        filedata = filedata.replace('JOB_NAME', name)
        bsub_file = os.path.join(os.getcwd(), '{}_oe_wbo_parent.lsf'.format(name))
        with open(bsub_file, 'w') as f:
            f.write(filedata)

        #submit job
        #stdin_file = open(bsub_file, 'r')
        #subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
        #stdin_file.close()

