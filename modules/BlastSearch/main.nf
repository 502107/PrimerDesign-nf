#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
=================================
= Primer Design Nextflow Script =
=================================
*/

params.baseDir = System.getProperty("user.dir")
params.hostRef = "${params.baseDir}/assets/refs/host/GCF_018294505.1.fna"
params.pgtRefsDir = "${params.baseDir}/assets/refs/pgt_all/*.fna"
params.pstRefsDir = "${params.baseDir}/assets/refs/pst_all/*.fna"
params.outputDir = "${params.baseDir}/results"
params.direction = "forward"

host = Channel.fromPath(params.hostRef)
pgt_refs = Channel.fromPath(params.pgtRefsDir)
pst_refs = Channel.fromPath(params.pstRefsDir)

process CombineRustReferences {
    publishDir "${params.outputDir}/combined_refs", mode: 'copy'
    
    input:
    path pgt_files
    path pst_files
    
    output:
    path 'pgt_all.fna', emit: pgt_all
    path 'pst_all.fna', emit: pst_all
    
    script:
    """
    cat ${pgt_files} > pgt_all.fna
    cat ${pst_files} > pst_all.fna
    """
}

rust_isol = Channel.from(
    tuple('pgt', '210', "${params.baseDir}/assets/refs/pgt210"),
    tuple('pst', '104', "${params.baseDir}/assets/refs/pst104")
)

include { BlastSearch as BlastSearch1 } from './modules/BlastSearch'
include { BlastSearch as BlastSearch2 } from './modules/BlastSearch'

workflow {
    combinedRefs = CombineRustReferences(pgt_refs.collect(), pst_refs.collect())
    pgt_all_ref = combinedRefs.pgt_all
    pst_all_ref = combinedRefs.pst_all
    
    genes_ch = ExtractValidGenes(rust_isol)

    pgt_genes_ch = genes_ch.filter { it[0] == 'pgt' } 
    pst_genes_ch = genes_ch.filter { it[0] == 'pst' } 

    pgt_blast_results_ch = BlastSearch1(pgt_genes_ch, pgt_refs) 
    pst_blast_results_ch = BlastSearch2(pst_genes_ch, pst_refs)

    blast_out_ch = pgt_blast_results_ch.mix(pst_blast_results_ch)
    
    // Ensure that BlastSearch1 and BlastSearch2 always emit their respective
    // (rust_type, path_to_blast_output) tuples. If no blast hits, path_to_blast_output can be an empty file.
    // This will ensure blast_out_ch.groupTuple() contains entries for both 'pgt' and 'pst'.
    conserved_ch = FindConservedSites(blast_out_ch.groupTuple())
    DesignPrimers(conserved_ch, host, pgt_all_ref, pst_all_ref)
}

process ExtractValidGenes {
    tag "${rust}"
    beforeScript "rm -f ${baseDir}/assets/refs/*/sorted_*.fna"
    executor 'slurm'
    cpus 1
    memory '4 GB'
    time '1h'

    publishDir "${params.outputDir}/genes", mode: 'copy'

    input:
    tuple val(rust), val(isol), val(isol_dir)

    output:
    tuple val(rust), path("${rust}_genes.fna")

    script:
    """
    rust_list="${params.baseDir}/${rust}_genes.txt"
    isol_gff="${isol_dir}/*.gff"
    isol_fa="${isol_dir}/*.fna"
    
    extract_genes_validate.py \$isol_fa \$isol_gff \$rust_list ${rust}_genes.fna
    """
}

process FindConservedSites {
    tag "${rust}"
    executor 'slurm'
    cpus 2
    memory '8 GB'
    time '2h'
    publishDir "${params.outputDir}/conserved_sites", mode: 'copy'

    input:
    tuple val(rust), path(blast_files) // blast_files will be a list of paths from groupTuple, staged by Nextflow

    output:
    tuple val(rust), path("${rust}_conserved_sites.fna")

    script:
    """
    source package 999eb878-6c39-444e-a291-e2e0a86660e6 # Load clustalw2 through the prokka package
    conserved_sites.py ${rust} "${params.outputDir}/blast_results" "${params.outputDir}/genes/${rust}_genes.fna" ${rust}_conserved_sites.fna
    """
}

process DesignPrimers {
    tag "${rust}"
    executor 'slurm'
    cpus 8
    memory '16 GB'
    time '4h'
    publishDir "${params.outputDir}/primers", mode: 'copy'

    input:
    tuple val(rust), path(conserved_sites)
    path host_ref
    path pgt_all_ref
    path pst_all_ref

    output:
    path "${rust}_primers.csv", optional: true

    script:
    def isol_dir = rust == "pgt" ? "${params.baseDir}/assets/refs/pgt210" : "${params.baseDir}/assets/refs/pst104"
    def actual_rust_ref = rust == "pgt" ? pgt_all_ref : pst_all_ref
    """
    source package 999eb878-6c39-444e-a291-e2e0a86660e6 # Load GNU Parallel through the prokka package

    mkdir -p tmp_${rust}_genes
    split_fasta.py ${conserved_sites} tmp_${rust}_genes/

    find tmp_${rust}_genes/ -type f -name "*.fna" | \
    parallel -j ${task.cpus} "bin/design_primers.py {} ${isol_dir} ${host_ref} ${params.direction} ${actual_rust_ref} {/.}_primers.tmp"

    if ls tmp_${rust}_genes/*_primers.tmp 1> /dev/null 2>&1; then
        cat tmp_${rust}_genes/*_primers.tmp > ${rust}_primers.csv
    else
        touch ${rust}_primers.csv
    fi
    """
}