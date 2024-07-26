process GENERATE_COMMANDS {

    tag "${meta.id}"
    container "google/deepvariant:1.5.0-gpu"
    label "process_single"

    input:
    each dv_args
    tuple val(meta),
          path(bam),
          path(bam_bai),
          path(ref),
          path(ref_fai)

    output:
    tuple val(meta),
          path(bam),
          path(bam_bai),
          path(ref),
          path(ref_fai),
          path("dv-make-examples.sh"),
          path("dv-call-variants-and-postprocess.sh")

    script:
    // declare optional flags if present in dv_args
    dv_opt_call_variants_extra_args = dv_args.call_variants_extra_args ? "--call_variants_extra_args ${dv_args.call_variants_extra_args}" : ""
    dv_opt_make_examples_extra_args = dv_args.make_examples_extra_args ? "--make_examples_extra_args ${dv_args.make_examples_extra_args}" : ""
    dv_opt_postprocess_variants_extra_args = dv_args.postprocess_variants_extra_args ? "--postprocess_variants_extra_args ${dv_args.postprocess_variants_extra_args}" : ""
    dv_opt_model_type = dv_args.model_type ? "--model_type ${dv_args.model_type}" : ""
    dv_opt_regions = dv_args.regions ? "--regions ${dv_args.regions}" : ""
    dv_opt_num_shards = dv_args.num_shards ? "--num_shards ${dv_args.num_shards}" : ""
    """
    # get commands and pipe to file
    run_deepvariant \\
        --dry_run \\
        --ref "${ref}" \\
        --reads "${bam}" \\
        --intermediate_results_dir . \\
        --output_vcf "${meta.id}.vcf.gz" \\
        ${dv_opt_model_type} ${dv_opt_regions} ${dv_opt_num_shards} \\
        ${dv_opt_call_variants_extra_args} ${dv_opt_make_examples_extra_args} ${dv_opt_postprocess_variants_extra_args} \\
        | grep time > cmds.sh

    # split commands
    # post process doesn't use GPU, but is inexpensive, so is lumped with call variants.
    head -n 1 cmds.sh > dv-make-examples.sh
    tail -n 2 cmds.sh > dv-call-variants-and-postprocess.sh
    chmod +x dv-make-examples.sh dv-call-variants-and-postprocess.sh
    """

}

process MAKE_EXAMPLES {

    tag "${meta.id}"
    container "google/deepvariant:1.5.0-gpu"
    cpus "${dv_args.num_shards ? dv_args.num_shards : 1}"
    memory "${dv_args.num_shards} GB" // doesn't needs *heaps* of memory
    label "process_long"

    input:
    each dv_args
    tuple val(meta),
          path(bam),
          path(bam_bai),
          path(ref),
          path(ref_fai),
          path(make_examples_script),
          path(call_variants_script)

    output:
    tuple val(meta),
          path(bam),
          path(bam_bai),
          path(ref),
          path(ref_fai),
          path(make_examples_script),
          path(call_variants_script),
          path("make_examples.*", type: "file")

    script:
    """
    . ${make_examples_script}
    """

}

process GPU_PART {

    tag "${meta.id}"
    container "google/deepvariant:1.5.0-gpu"
    publishDir "${params.outdir}", pattern: "*.html", mode: "copy" // only save the report
    label "process_gpu"

    input:
    tuple val(meta),
          path(bam),
          path(bam_bai),
          path(ref),
          path(ref_fai),
          path(make_examples_script),
          path(call_variants_script),
          path(dv_examples)

    output:
    tuple val(meta),
          path(bam),
          path(bam_bai),
          path(ref),
          path(ref_fai),
          path("${meta.id}.vcf.gz"),
          path("${meta.id}.vcf.gz.tbi"),
          path("${meta.id}.visual_report.html")

    script:
    """
    . ${call_variants_script}
    """
}

workflow separated_deepvariant {

    take:
        deepvariant_params
        reads
        ref
    main:
        // combine channels and feed into deepvariant
        reads
            .combine(ref.map{tuple(it[1], it[2])})
            .set {ch_dv_input}

        ch_cmds = GENERATE_COMMANDS(deepvariant_params, ch_dv_input)
        MAKE_EXAMPLES(deepvariant_params, ch_cmds)
            | GPU_PART
            | set{results}
    emit:
        results
}

// set params to a falsey default value
params.model_type = ''
params.regions = ''
params.num_shards = ''
params.make_examples_extra_args = ''
params.call_variants_extra_args = ''
params.postprocess_variants_extra_args = ''
params.bams_dir = ''
params.ref = ''

def helpMessage() {
    log.info"""
    Wrapper pipeline around Google GPU DeepVariant run_deepvariant script.

    Usage:
        nextflow run main.nf --bams_dir <path> --ref <path> [options]

    Required Arguments:
        --bams_dir: Aligned, sorted, indexed BAM file containing the reads we
          want to call. Should be aligned to a reference genome compatible with --ref.
        --ref: Genome reference to use. Must have an associated FAI index as
          well. Supports text or gzipped references. Should match the reference used
          to align the BAM file provided to --reads.

    Optional Arguments:
        --call_variants_extra_args: A comma-separated list of flag_name=flag_value.
          "flag_name" has to be valid flags for call_variants.py. If the flag_value is
          boolean, it has to be flag_name=true or flag_name=false.
        --make_examples_extra_args: A comma-separated list of flag_name=flag_value.
          "flag_name" has to be valid flags for make_examples.py. If the flag_value is
          boolean, it has to be flag_name=true or flag_name=false.
        --num_shards: Number of shards for make_examples step.
          (default: '1')
          (an integer)
        --postprocess_variants_extra_args: A comma-separated list of
          flag_name=flag_value. "flag_name" has to be valid flags for
          postprocess_variants.py. If the flag_value is boolean, it has to be
          flag_name=true or flag_name=false.
        --regions: Space-separated list of regions we want to process.
          Elements can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE
          files.""".stripIndent()

}

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }
    ch_dv_args = channel.from(
        [
            model_type: params.model_type,
            regions: params.regions,
            num_shards: params.num_shards,
            make_examples_extra_args: params.make_examples_extra_args,
            call_variants_extra_args: params.call_variants_extra_args,
            postprocess_variants_extra_args: params.postprocess_variants_extra_args
        ]
    )
    ch_samples = channel.fromPath("${params.bams_dir}/*.bam")
        .map{tuple([id: it.baseName], it, "${it}.bai")}
    ch_ref = channel.fromPath("${params.ref}", checkIfExists: true)
        .map{tuple([id: it.baseName], it, "${it}.fai")}
    separated_deepvariant(ch_dv_args, ch_samples, ch_ref)
}
