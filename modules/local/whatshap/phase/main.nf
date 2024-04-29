process WHATSHAP_PHASE {

    tag "$meta.id"
    label "process_single"

    // 1.4 for compatibility
    conda "bioconda::whatshap=1.4"
    container "quay.io/biocontainers/whatshap:1.4--py39hc16433a_0"

    input:
    tuple val(meta), path(bam), path(bam_idx)
    tuple val(meta_vcf), path(vcf_gz_PASS), path(vcf_gz_PASS_idx)
    tuple val(meta_ref), path(reference), path(reference_idx)

    output:
    tuple val(meta), path("*.phased.vcf.gz"), path("*.phased.vcf.gz.tbi")

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