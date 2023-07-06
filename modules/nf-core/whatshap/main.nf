process WHATSHAP {
    tag "haplotagging.."
    label 'regular'
    cpus 8
    memory '32 GB'
    time '48h'
    queue 'regular'
    executor 'slurm'
    publishDir "$params.outdir/whatshap", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::whatshap=1.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:1.7--py39hc16433a_0' :
        'quay.io/biocontainers/whatshap:1.7--py39hc16433a_0' }"
    
   

    input:
    tuple val(lib), path(deepvariant_vcf)
    tuple val(lib), path(deepvariant_idx)
    tuple val(lib), path(ontfile_bam)
    tuple val(lib), path(ontfile_idx)
    
    
    output:
    tuple val(lib), path("*.bam"), emit: bam

    script:
    """

    whatshap haplotag \\
        --ignore-read-groups \\
        --output ./${lib}.hp.bam \\
        --reference $params.fasta \\
        $deepvariant_vcf \\
        $ontfile_bam \
       
    """





}
