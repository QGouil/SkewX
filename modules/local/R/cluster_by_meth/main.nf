process R_CLUSTERBYMETH {

    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}/cluster_by_meth", mode: "copy", pattern: "*_CGIX_skew.tsv.gz"
    publishDir "${params.outdir}/cluster_by_meth", mode: "copy", pattern: "*_CGIX_clustered_reads.tsv.gz"

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://qgouil/skewx/skewx-r:0.2' :
        'ghcr.io/qgouil/skewx-r:0.2' }"

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
