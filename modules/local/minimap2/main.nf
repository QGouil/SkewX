process MINIMAP2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::minimap2=2.17" // find environment or container that has both samtools and minimap2

    input:
    tuple val(meta), path(input_file) // meta keeps track of individual + their tissue samples
    tuple val(ref_basename), path(ref_fasta), path(ref_faidx) // ch_reference

    output:
    tuple val(meta), path("*mapped.bam")

    script:
    """
    #convert to fastq, keeping modbam tags
    samtools fastq -@ 12 -T '*' -o ${meta.id}.fq.gz ${input_file}
    #if ont reads
    if [ "$params.lrs" == "ont" ]; then
        minimap2 -ax map-ont -t 12 ${ref_fasta} ${meta.id}.fq.gz | samtools sort > ${meta.id}_mapped.bam
    elif [ "$params.lrs" == "pacbio" ]; then
        minimap2 -ax map-pb -t 12 ${ref_fasta} ${meta.id}.fq.gz | samtools sort > ${meta.id}_mapped.bam
    fi
    """
}