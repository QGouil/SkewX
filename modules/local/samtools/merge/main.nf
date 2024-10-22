process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(input_files), path(bais) // meta keeps track of individual + their tissue samples

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    """
    # make new merged bam
    samtools \\
        merge \\
        --threads ${task.cpus} \\
        ${meta.id}_tmp.bam \\
        ${input_files}
    
    # remove read groups
    samtools view -H ${meta.id}_tmp.bam | \\
        grep -v "^@RG" | \\
        samtools reheader - ${meta.id}_tmp.bam > \\
        ${meta.id}.bam
    """
}
