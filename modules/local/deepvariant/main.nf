process DEEPVARIANT {
    tag "$meta.id"
    label 'process_gpu'
    publishDir "${params.outdir}", pattern: "*.html", mode: "copy" // only save the report
    container "google/deepvariant:1.5.0-gpu"

    input:
    each region
    each model
    tuple val(meta), path(input_reads), path(input_reads_idx) // meta keeps track of individual + their tissue samples
    tuple val(meta_ref), path(reference), path(reference_idx)

    output:
    tuple val(meta), path(input_reads), path(input_reads_idx), path("*.vcf.gz") // keep mutations attached to BAM
    tuple val(meta), path("*.html"), emit: html

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type ${model} \\
        --ref ${reference} \\
        --reads ${input_reads} \\
        --regions ${region} \\
        --output_vcf ${meta.id}.vcf.gz \\
        --num_shards ${task.cpus}
    """
}