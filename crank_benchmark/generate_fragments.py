import fragmenter
import os

# Generate fragments. This will generate a dictionary that maps fragments to parent molecule. It can be written out as a
# JSON file.
fragments = fragmenter.generate_fragments('input_molecules.smi', remove_map=False, json_filename='fragments.json')

# Find rotatable bonds in fragments and generate specs for crank jobs. This step removes duplicate fragments.
crank_jobs = fragmenter.fragment_to_torsion_scan(fragments, json_filename='crank_jobs.json')

# Launch crank jobs
for fragment in crank_jobs:
    # Generate input files for crank-launch
    for job in crank_jobs[fragment]['crank_torsion_drives']:
        path, inputfile, dihedraltxt = fragmenter.to_crank_input(crank_jobs[fragment], crank_job=job)
        # Make necessary changes to submit script and submit crank job
        with open('submit_dummy', 'r') as f:
            filedata = f.read()
        JOB_NAME = '{}'.format(crank_jobs[fragment]['canonical_isomeric_SMILES'])
        filedata = filedata.replace('JOB_NAME', '{}'.format(JOB_NAME))
        filedata = filedata.replace('INPUTFILE', '{}'.format(inputfile))
        filedata = filedata.replace('DIHEDRALS', '{}'.format(dihedraltxt))
        outfile = inputfile.replace('dat', 'out')
        filedata = filedata.replace('OUTFILE', '{}'.format(outfile))

        bsub_file = os.path.join(path, 'submit.sh')
        with open(bsub_file, 'w') as f:
            f.write(filedata)
        #change directory and submit job
        os.chdir(path)
        os.system('bsub < {}'.format(bsub_file))
        os.chdir('..')


