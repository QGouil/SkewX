process WHATSHAP_PHASE {

    tag "$meta.id"
    label "process_low"

    conda "bioconda::whatshap=2.2"
    container "quay.io/biocontainers/whatshap:2.2--py310h0dbaff4_1"

    input:
    tuple val(meta), path(bam)
    tuple val(meta_idx), path(bam_idx)
    tuple val(meta_vcf), path(vcf_gz_PASS)
    tuple val(meta_vcf_idx), path(vcf_gz_PASS_idx)
    tuple val(meta_ref), path(reference)
    tuple val(meta_ref_idx), path(reference_idx)

    output:
    tuple val(meta), path("*.phased.vcf.gz"), emit: gz
    tuple val(meta), path("*.phased.vcf.gz.tbi"), emit: tbi

    script:
    """
    whatshap phase \\
        --ignore-read-groups \\
        --output "${meta.id}.phased.vcf.gz" \\
        --reference "${reference}" \\
        "${vcf_gz_PASS}" \\
        "${bam}"
    
    # index the compressed vcf
    tabix -p vcf "${meta.id}.phased.vcf.gz"
    """
}