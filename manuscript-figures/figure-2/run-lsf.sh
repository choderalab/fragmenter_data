#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00
#
# Specify node group
#BSUB -q cpuqueue
#
#
# nodes: number of nodes and GPU request
#BSUB -n 1
#BSUB -R "rusage[mem=2] span[hosts=1]"
#
# job name (default = name of script file)
#BSUB -J "drug-bank-wbo[1-10]"

# Make sure to run bashrc
source $HOME/.bashrc

export NJOBS=10

# Activate conda env
conda activate qcf

# Process dataset fragment
python generate_wbo.py $LSB_JOBINDEX $NJOBS drugbank-wbo-${LSB_JOBINDEX}.oeb