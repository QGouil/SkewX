process MOSDEPTH {

    tag "$meta.id"
    label "process_low"

    conda "bioconda::mosdepth=0.3.6"
    container 'quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0'

    input:
    tuple val(meta), path(haplotagged_bam), path(sample_bam_idx)

    output:
    tuple val(meta), path('*.bed.gz'), path('*.bed.gz.csi')
    tuple val(meta), path('*mosdepth.global.dist.txt'), path('*mosdepth.region.dist.txt'), path('*mosdepth.summary.txt')

    script:
    """
    mosdepth \\
        -t ${task.cpus} \\
        -n \\
        --fast-mode \\
        --by 1000000 \\
        "${meta.id}_${meta.sample}" \\
        "${haplotagged_bam}"
    """
}