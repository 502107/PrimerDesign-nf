process DesignPrimers {
    tag "${rust}"
    executor 'slurm'
    cpus 32
    memory '200 GB'
    publishDir "${params.outputDir}/primers", mode: 'copy'

    input:
    tuple val(rust), path(conserved_sites), val(padding)
    path host_ref
    path rust_all
    path blast_out

    output:
    path "${rust}_primers.csv", optional: true

    script:
    def isol_dir = rust == "pst" ? "${params.baseDir}/assets/refs/pst130" : "${params.baseDir}/assets/refs/pgt210"
    """
    source package 999eb878-6c39-444e-a291-e2e0a86660e6 # Load GNU Parallel through the prokka package

    mkdir -p tmp_${rust}_genes
    split_fasta.py ${conserved_sites} tmp_${rust}_genes/

    find tmp_${rust}_genes/ -type f -name "*.fna" | \
    parallel -j ${task.cpus} "${params.baseDir}/bin/design_primers.py {} ${rust} ${isol_dir} ${host_ref} ${blast_out} ${rust_all} {}_primers.tmp ${padding} ${params.direction}"

    if ls tmp_${rust}_genes/*_primers.tmp 1> /dev/null 2>&1; then
        (head -qn 1 tmp_${rust}_genes/*_primers.tmp | uniq && tail -qn 1 tmp_${rust}_genes/*_primers.tmp | sort -rnk 2) > ${rust}_primers.csv
    else
        touch ${rust}_primers.csv
    fi
    """
}