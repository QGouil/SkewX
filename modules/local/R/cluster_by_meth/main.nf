process R_CLUSTERBYMETH {
    
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/../environment.yml"

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