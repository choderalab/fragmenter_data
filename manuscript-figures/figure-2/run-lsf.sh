#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 4:00
#
# Specify node group
#BSUB -q cpuqueue
#
#
# nodes: number of nodes and GPU request
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"
#
# job name (default = name of script file)
#BSUB -J "drug-bank-wbo[1-21]"
#BSUB -o /home/chayas/job_output/drugbank-wbo.stdout
#BSUB -eo /home/chayas/job_output/drugbank-wbo.stderr

# Make sure to run bashrc
source $HOME/.bashrc

export NJOBS=20

# Activate conda env
conda activate qcf-manager-openff

# Process dataset fragment
python generate_wbo.py $LSB_JOBINDEX $NJOBS drugbank-wbo-${LSB_JOBINDEX}.oeb