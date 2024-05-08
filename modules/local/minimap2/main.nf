process MINIMAP2 {
    tag "$meta.id"
    label 'process_high'
    //conda "bioconda::minimap2=2.17" // find environment or container that has both samtools and minimap2

    input:
    tuple val(meta), path(input_file) // meta keeps track of individual + their tissue samples
    tuple val(ref_basename), path(ref_fasta), path(ref_faidx) // ch_reference

    output:
    tuple val(meta), path("${meta.id}_${meta.sample}_mapped.bam")

    script:
    // set the minimap2 preset appropriately
    if (params.lrs == "ont") {
        minimap2_preset = "map-ont"
    } else if (params.lrs == "pacbio") {
        minimap2_preset = "map-pb"
    } else {
        error "--lrs must be \"ont\" or \"pacbio\"!"
    }
    """
    # combine fastq/minimap2/sort to overlap IO with computations
    samtools fastq -T 'MM,ML' ${input_file} | \\
        minimap2 -ax $minimap2_preset -t $task.cpus "${ref_fasta}" - | \\
        samtools sort -@ $task.cpus > "${meta.id}_${meta.sample}_mapped.bam"
    """
}
