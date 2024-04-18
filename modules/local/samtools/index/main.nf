process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    container "${params.deepvariant_container}"

    input:
    tuple val(meta), path(input_file) // meta keeps track of individual + their tissue samples

    output:
    tuple val(meta), path("*.bam.bai"), emit: bai

    script:
    """
    # index merged bams
    samtools index -@ ${task.cpus} ${input_file}
    """
}