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
params.padding = "500"

host = Channel.fromPath(params.hostRef)
pgt_refs = Channel.fromPath(params.pgtRefsDir)
pst_refs = Channel.fromPath(params.pstRefsDir)

rust_isol = Channel.from(
    tuple('pgt', '210', "${params.baseDir}/assets/refs/pgt210", params.padding),
    tuple('pst', '130', "${params.baseDir}/assets/refs/pst130", params.padding)
)

include { BlastSearch as BlastSearchPgt } from './modules/BlastSearch'
include { BlastSearch as BlastSearchPst } from './modules/BlastSearch'
include { DesignPrimers as DesignPrimersPgt } from './modules/DesignPrimers'
include { DesignPrimers as DesignPrimersPst } from './modules/DesignPrimers'

workflow {
    pgt_all = pgt_refs.collectFile(name: 'pgt_all.fna', newLine: true)
    pst_all = pst_refs.collectFile(name: 'pst_all.fna', newLine: true)
    
    genes_ch = ExtractValidGenes(rust_isol)

    pgt_genes_ch = genes_ch.filter { it[0] == 'pgt' } 
    pst_genes_ch = genes_ch.filter { it[0] == 'pst' } 

    pgt_blast_results_ch = BlastSearchPgt(pgt_genes_ch, pgt_refs) 
    pst_blast_results_ch = BlastSearchPst(pst_genes_ch, pst_refs)

    blast_out_ch = pgt_blast_results_ch.mix(pst_blast_results_ch)
    
    conserved_ch = FindConservedSites(blast_out_ch.groupTuple())
    // Split conserved_ch by rust type
    pgt_conserved_ch = conserved_ch.filter { it[0] == 'pgt' }
    pst_conserved_ch = conserved_ch.filter { it[0] == 'pst' }

    def blast_output_dir = file("${params.outputDir}/blast_results")
    
    DesignPrimersPgt(pgt_conserved_ch.map{ rust_val, sites_path -> tuple(rust_val, sites_path, params.padding) }, host, pgt_all, blast_output_dir)
    DesignPrimersPst(pst_conserved_ch.map{ rust_val, sites_path -> tuple(rust_val, sites_path, params.padding) }, host, pst_all, blast_output_dir)

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
    tuple val(rust), val(isol), val(isol_dir), val(padding)

    output:
    tuple val(rust), path("${rust}_genes.fna")

    script:
    """
    rust_list="${params.baseDir}/${rust}_genes.txt"
    isol_gff="${isol_dir}/*.gff"
    isol_fa="${isol_dir}/*.fna"
    
    extract_genes_validate.py \$isol_fa \$isol_gff \$rust_list ${rust}_genes.fna ${padding}
    """
}

process FindConservedSites {
    tag "${rust}"
    executor 'slurm'
    cpus 16
    memory '50 GB'
    time '20h'
    publishDir "${params.outputDir}/conserved_sites", mode: 'copy'

    input:
    tuple val(rust), path(blast_files)

    output:
    tuple val(rust), path("${rust}_conserved_sites.fna")

    script:
    """
    source package 999eb878-6c39-444e-a291-e2e0a86660e6 # Load clustalw2 through the prokka package
    conserved_sites.py ${rust} "${params.outputDir}/blast_results" "${params.outputDir}/genes/${rust}_genes.fna" ${rust}_conserved_sites.fna
    """
}