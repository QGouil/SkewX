process FILTER_PASS {
    tag "$meta.id"
    label "process_low"
    // publishDir
    conda "bioconda::bcftools=1.17"
    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"

    input:
    tuple val(meta), path(bam), path(bam_idx), path(mutations)

    output:
    tuple val(meta), path(bam), path(bam_idx), path("*_PASS.vcf.gz"), path("*_PASS.vcf.gz.tbi")

    script:
    """
    bcftools view -f PASS "${mutations}" > "${meta.id}_PASS.vcf"
    bgzip "${meta.id}_PASS.vcf"
    tabix -p vcf "${meta.id}_PASS.vcf.gz"
    """
}