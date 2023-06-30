process NANOCOMP {
    tag "QC"
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
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png") , emit: png
    tuple val(meta), path("*.txt") , emit: txt
    tuple val(meta), path("*.log") , emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_file = ("$ontfile".endsWith(".bam")) ? "--bam ${ontfile}" : ''
    """
    NanoComp \\
        $args \\
        -t 4 \\
        -o ./ \\
        $input_file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocomp: \$(echo) \$(NanoComp --version 2>&1) | sed 's/^.*NanoComp //; s/ .*//')
    END_VERSIONS
    """
}



