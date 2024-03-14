#!/bin/bash -e
#SBATCH --mem 20G
#SBATCH -J mainPD
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH -p jic-short,nbi-short,jic-long,nbi-long
#SBATCH -o ./mainPD_%j.out
#SBATCH -e ./mainPD_%j.err

source ~/mambaforge/etc/profile.d/conda.sh
source ~/mambaforge/etc/profile.d/mamba.sh
mamba activate nextflow

./main_splt.nf -resume

