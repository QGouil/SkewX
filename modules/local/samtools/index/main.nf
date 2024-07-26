process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/whatshap", mode: 'copy', pattern: '*hp.bam.bai'

    conda "bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(input_file) // meta keeps track of individual + their tissue samples

    output:
    tuple val(meta), path("*.bam", includeInputs: true), path("*.bam.bai")

    script:
    """
    # index merged bams
    samtools index -@ ${task.cpus} ${input_file}
    """
}
