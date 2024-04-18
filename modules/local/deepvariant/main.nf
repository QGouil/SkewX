process DEEPVARIANT {
    tag "$meta.id"
    label 'process_gpu'
    publishDir "${params.outdir}", pattern: "*.html", mode: "copy" // only save the report
    container "google/deepvariant:1.5.0-gpu"

    input:
    val(model)
    tuple val(meta), path(input_reads) // meta keeps track of individual + their tissue samples
    tuple val(meta_idx), path(input_reads_index)
    tuple val(meta_ref), path(reference)
    tuple val(meta_ref_idx), path(reference_index)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: gz
    tuple val(meta), path("*.html"), emit: html

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type ${model} \\
        --ref ${reference} \\
        --reads ${input_reads} \\
        --regions chrX \\
        --output_vcf ${meta.id}.vcf.gz \\
        --num_shards ${task.cpus}
    """
}