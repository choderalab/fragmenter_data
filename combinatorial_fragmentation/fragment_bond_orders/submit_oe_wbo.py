import os
import subprocess
jobs_to_run = ['Larotrectinib', 'Lenvatinib', 'Nilotinib', 'Palbociclib', 'Pazopanib', 'Ponatinib', 'Ribociclib',
               'Sunitinib', 'Tofacitinib', 'Vemurafenib']

for name in jobs_to_run:
    with open('oe_wbo.lsf', 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('JOB_NAME', name)
    bsub_file = os.path.join(os.getcwd(), '{}_oe_wbo.lsf'.format(name))
    with open(bsub_file, 'w') as f:
        f.write(filedata)

    #submit job
    stdin_file = open(bsub_file, 'r')
    subprocess.Popen(['bsub'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin_file)
    stdin_file.close()

