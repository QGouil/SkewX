process FILTER_PASS {
    tag "$meta.id"
    label "process_low"
    // publishDir
    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"

    input:
    tuple val(meta), path(mutations)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: gz
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi

    script:
    """
    bcftools view -f PASS "${mutations}" > "${meta.id}_PASS.vcf"
    bgzip "${meta.id}_PASS.vcf"
    tabix -p vcf "${meta.id}_PASS.vcf.gz"
    """
}