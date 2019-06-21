import os
import subprocess
import glob
# jobs_to_run = ['Larotrectinib', 'Lenvatinib', 'Nilotinib', 'Palbociclib', 'Pazopanib', 'Ponatinib', 'Ribociclib',
#                'Sunitinib', 'Tofacitinib', 'Vemurafenib']

jobs_to_run = glob.glob('../fragment/validation_set/fragments*')

for file in jobs_to_run:
    name = file.split('/')[-1].split('.')[0].split('_')
    name = name[1] + '_' +  name[2]
    print(name)
    with open('oe_wbo.lsf', 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('JOB_NAME', name)
    bsub_file = os.path.join(os.getcwd(), 'validation_set/{}_oe_wbo.lsf'.format(name))
    with open(bsub_file, 'w') as f:
        f.write(filedata)

    #submit job
    os.chdir(os.path.join(os.getcwd(), 'validation_set'))
    stdin_file = open(bsub_file, 'r')
    subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
    stdin_file.close()
    os.chdir('..')

