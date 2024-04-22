process MOSDEPTH_PLOTDIST {

    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}", mode: "copy"

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(plot_dist_script)
    tuple val(meta), path(global_dist)

    output:
    tuple val(meta), path("${meta.id}.dist.html"), emit: html

    script:
    """
    python3 "${plot_dist_script}" *.global.dist.txt
    mv dist.html ${meta.id}.dist.html
    """
    
}