process BlastSearch {
    tag "${rust_of_gene}_vs_${individual_ref_file.getBaseName()}"
    executor 'slurm'
    cpus 4
    memory '8 GB'
    time '2h'
    publishDir "${params.outputDir}/blast_results", mode: 'copy'

    input:
    tuple val(rust_of_gene), path(genes_file)
    each individual_ref_file

    output:
    tuple val(rust_of_gene), path("${rust_of_gene}_vs_${individual_ref_file.getBaseName()}_blast.out")

    script:
    """
    blastn -query ${genes_file} -subject ${individual_ref_file} \\
           -outfmt "6 qseqid sseqid pident bitscore sstart send sseq" \\
           -out ${rust_of_gene}_vs_${individual_ref_file.getBaseName()}_blast.out

    touch ${rust_of_gene}_vs_${individual_ref_file.getBaseName()}_blast.out
    """
}
