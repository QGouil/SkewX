process NANOCOMP {
    tag "performing QC.."
    label 'regular'
    cpus 8
    memory '32 GB'
    time '48h'
    queue 'regular'
    executor 'slurm'

    conda (params.enable_conda ? 'bioconda::nanocomp=1.22.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp:1.22.0--pyhdfd78af_0' :
        'quay.io/biocontainers/nanocomp:1.22.0--pyhdfd78af_0' }"
    
    publishDir "$params.outdir/nanocomp", mode: 'copy'

    input:
    tuple val(lib), path(ontfile_bam)

    output:
    tuple val(lib), path("*.html"), emit: html
    tuple val(lib), path("*.png") , emit: png
    tuple val(lib), path("*.txt") , emit: txt
    tuple val(lib), path("*.log") , emit: log

    script:
    """
    NanoComp \\
        -o ./ \\
	--verbose \\
        --bam $ontfile_bam \
    """
}



