process R_CLUSTERBYMETH {
    
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/../environment.yml"
    container "oras://ghcr.io/wehi-researchcomputing/skewx-r:0.1"

    input:
    tuple val(meta), path(bam), path(bam_idx), path(hpreads)
    tuple val(meta_bed), path(cgix_bed)

    output:
    tuple val(meta), path("${meta.id}_${meta.sample}_CGIX_clustered_reads.tsv.gz"), path("${meta.id}_${meta.sample}_CGIX_skew.tsv.gz")

    script:
    """
    cluster_by_meth.R "${meta.id}_${meta.sample}" "${bam}" "${cgix_bed}" "${hpreads}" ${task.cpus}
    """

}