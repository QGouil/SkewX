process WHATSHAP_HAPLOTAG {

    tag "$meta.id"
    label "process_low"

    conda "bioconda::whatshap=2.2"
    container "quay.io/biocontainers/whatshap:2.2--py310h0dbaff4_1"

    input:
    tuple val(meta), path(bam), path(bam_idx)
    tuple val(meta_vcf), path(vcf_gz), path(vcf_gz_idx)
    tuple val(meta_ref), path(ref), path(ref_idx)

    output:
    tuple val(meta), path('*.hp.bam')

    script:
    """
    whatshap haplotag \\
        --ignore-read-groups \\
        --output "${meta.id}.hp.bam" \\
        --reference "${ref}" \\
        "${vcf_gz}" \\
        "${bam}"
    """
}