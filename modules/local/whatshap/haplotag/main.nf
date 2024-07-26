process WHATSHAP_HAPLOTAG {

    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}/whatshap", mode: "copy", pattern: "*.hp.bam"

    // 1.4 for compatibility
    conda "bioconda::whatshap=1.4"
    container "quay.io/biocontainers/whatshap:1.4--py39hc16433a_0"

    input:
    tuple val(meta), path(bam), path(bam_idx), path(vcf_gz), path(vcf_gz_idx)
    tuple val(meta_ref), path(ref), path(ref_idx)

    output:
    tuple val(meta), path('*.hp.bam')

    script:
    """
    whatshap haplotag \\
        --ignore-read-groups \\
        --output "${bam.baseName}.hp.bam" \\
        --reference "${ref}" \\
        "${vcf_gz}" \\
        "${bam}"
    """
}
