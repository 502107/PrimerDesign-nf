#!/bin/bash -e
#SBATCH --mem 20G
#SBATCH -J PDnf
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH -p jic-long,nbi-long
#SBATCH -o ./PD_%j.out
#SBATCH -e ./PD_%j.err
#SBATCH -t 31-00:00

source ~/mambaforge/etc/profile.d/conda.sh
source ~/mambaforge/etc/profile.d/mamba.sh
mamba activate nextflow

rust=$1
ref_isol=$2
fasta=$3
host=$4
blast_results=$5
all_genes=$6
direction=$7
out=$8
basedir=$9

echo "Submitted $out"

python $basedir/bin/design_primers_splt.py $rust $ref_isol $fasta $host $blast_results $all_genes $out $direction
