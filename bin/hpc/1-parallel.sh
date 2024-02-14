#!/bin/bash
#SBATCH --mem 1G
#SBATCH -J clerkPD
#SBATCH -c 1
#SBATCH -p jic-long,nbi-long,RG-Diane-Saunders
#SBATCH -o clerk_%j.out
#SBATCH -e clerk_%j.err

rust=$1
ref_isol=$2
host=$3
genes=$4
blast_results=$5
all_refs=$6
out_dir=$7
basedir=$8

mkdir -p $out_dir
# Reading two lines at a time
ls $basedir/$rust-tmp | cut -d'.' -f1 | while read ID 
do
    fasta=$basedir/$rust-tmp/$ID.fna

    Jobs=$(squeue -u $USER | grep 'PDnf' | wc -l)
    while [ $Jobs -gt 50 ]
    do
        sleep 5
        printf "."
        Jobs=$(squeue -u $USER | grep "PDnf" | wc -l)
    done
    echo "$ID submitted as sbatch $basedir/bin/hpc/2-submit.slurm $rust $ref_isol $fasta $host $blast_results $all_refs $genes $out_dir/${ID}_primer.csv"
    sbatch $basedir/bin/hpc/2-submit.slurm $rust $ref_isol $fasta $host $blast_results $all_refs $genes $out_dir/${ID}_primer.csv $basedir
done

