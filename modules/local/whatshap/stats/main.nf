process WHATSHAP_STATS {

    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}", mode: "copy"
    
    // 1.4 for compatibility
    conda "bioconda::whatshap=1.4"
    container "quay.io/biocontainers/whatshap:1.4--py39hc16433a_0"

    input:
    tuple val(meta), path(vcf_gz), path(vcf_gz_idx)

    output:
    tuple val(meta), path("${meta.id}_whatshapstats.txt"), path("${meta.id}_blocks.tsv")

    script:
    """
    whatshap stats --block-list "${meta.id}_blocks.tsv" "${vcf_gz}" > "${meta.id}_whatshapstats.txt"
    """
}