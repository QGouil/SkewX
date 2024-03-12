process SNIFFLES2 {
    tag "sniffing for SVs..."
    label "sliffles2"
    cpus 8
    memory '32 GB'
    time '48h'
    queue 'regular'
    executor 'slurm'
    publishDir "$params.outdir/sniffles2", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::sniffles=2.0.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.0.7--pyhdfd78af_0' :
        'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0' }"



    input:
    tuple val(lib), path(ontfile_bam)
    tuple val(lib), path(ontfile_idx)
    

    output:
        path "*.sniffles.vcf", emit: vcf
    
    script:
        def sniffles_args = params.sniffles_args ?: ''
    """
    sniffles \
        --threads 8 \
        --sample-id ${lib} \
        --output-rnames \
        --cluster-merge-pos $params.cluster_merge_pos \
        --input $ontfile_bam \
        --tandem-repeats ${params.tr_bed} \
        $sniffles_args \
        --vcf ${lib}.sniffles.vcf
    sed '/.:0:0:0:NULL/d' ${lib}.sniffles.vcf > tmp.vcf
    mv tmp.vcf ${lib}.sniffles.vcf
    """
}
