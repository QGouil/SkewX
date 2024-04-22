process SAMTOOLS_VIEWHP {
    
    tag "$meta.id"
    label "process_low"

    conda "bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(bam), path(bam_idx)
    tuple val(meta_bed), path(bed)

    output:
    tuple val(meta), path("${bam.baseName}_CGIX_hpreads.tsv.gz")

    shell:
    '''
    samtools view -d HP -L "!{bed}" "!{bam}" \\
        | awk -F '\\t' 'BEGIN{OFS="\\t"} {for (i=12; i<=NF; ++i) { if ($i ~ "^HP:i:|^PS:i:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; print $1, td["HP"], td["PS"]}' \\
        | gzip \\
        > "!{bam.baseName}_CGIX_hpreads.tsv.gz"
    '''

}