process NANOCOMP {

    tag "$meta.id"
    label "process_low"

    conda "bioconda::nanocomp=1.23.1"
    container "quay.io/biocontainers/nanocomp:1.23.1--pyhdfd78af_0"

    input:
    tuple val(meta), path(bams), path(bams_idx)

    output:
    tuple val(meta), path("*.html")

    script:
    """
    NanoComp \\
        -t ${task.cpus} \\
        --bam "${bams.join("\" \"")}" \\
        --prefix "${meta.id}_" \\
    """

}