#!/bin/bash -e
#SBATCH --mem 20G
#SBATCH -J mainPD
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH -p jic-long,nbi-long,RG-Diane-Saunders
#SBATCH -o ./mainPD_%j.out
#SBATCH -e ./mainPD_%j.err
#SBATCH -t 31-00:00

source ~/mambaforge/etc/profile.d/conda.sh
source ~/mambaforge/etc/profile.d/mamba.sh
mamba activate nextflow

./main.nf -resume

