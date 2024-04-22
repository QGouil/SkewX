#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rrms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rrms

    Website: https://nf-co.re/rrms
    Slack  : https://nfcore.slack.com/channels/rrms
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nanopore reads PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
if (params.input) { 
    ch_input = Channel.fromPath(params.input, checkIfExists: true) 
} else { exit 1, 'Input sample sheet not specified!' }
if (params.reference) { 
    ch_reference = Channel.fromPath("${params.reference}", checkIfExists: true).map{
        it -> tuple(id: it.baseName, it)
    }
    ch_reference_index = Channel.fromPath("${params.reference}.fai", checkIfExists: true).map{
        it -> tuple(id: it.baseName, it)
    }
} else { exit 1, "Reference FASTA not specified!" }
if (params.cgi_bedfile) {
    ch_cgibed = Channel.fromPath(params.cgi_bedfile, checkIfExists: true)
} else { exit 1, "CGI Bed file not specified!"}
if (!(params.deepvariant_model in ["WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA"])) {
    exit 1, "DeepVariant model must be one of WGS, WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA"
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { RRMS } from './workflows/rrms'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define processes and modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include {INPUT_CHECK} from './subworkflows/local/input_check.nf'
include {SAMTOOLS_MERGE} from "./modules/local/samtools/merge/main.nf"
// two processes needed as two sets of bams are indexed.
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_SAMPLES} from "./modules/local/samtools/index/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED} from "./modules/local/samtools/index/main.nf"
include {DEEPVARIANT} from "./modules/local/deepvariant/main.nf"
/*
process dorado_mod_basecall {
    debug true
    label 'gpu'
    cpus 16
    memory '120GB'
    time '48h'
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A100:1 --cpus-per-task=16 --qos=bonus'
    publishDir "$params.outdir/sup_5mCG_5hmCG_alignments", mode: 'copy'
    module 'samtools/1.19.2'
    module 'dorado/0.5.2'


    input:
        tuple val(lib), path(raw_dir)
    output:
        tuple val(lib), path("*.bam"), emit: bam
        tuple val(lib), path("*.bam.bai"), emit: bai
    script:
        """
        dorado basecaller "$params.dorado_config" ${raw_dir} --reference "$params.fasta" --modified-bases 5mCG_5hmCG | samtools sort -@ 8 > ${lib}_sup_5mCG_5hmCG.CHM13v2.bam
        samtools index -@ 8 ${lib}_sup_5mCG_5hmCG.CHM13v2.bam

        """
} */

process minimap2 {
    tag "Aligning reads.."
    label 'align'
    memory '80 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    clusterOptions '--qos=bonus'
    publishDir "$params.outdir/minimap2", mode: 'copy'
    module 'minimap2/2.17'
    module 'samtools/1.19.2'

    input:
        tuple val(lib), path(ontfile)
    output:
        tuple val(lib), path("*.bam"), emit: bam
        tuple val(lib), path("*.bam.bai"), emit: bai
    script:
        """
        minimap2 -ax map-ont -t 8 "$params.fasta" $ontfile | samtools sort -@ 8 > ${lib}_sup_5mCG_5hmCG.CHM13v2.bam
        samtools index -@ 8 ${lib}_sup_5mCG_5hmCG.CHM13v2.bam
        """

}

process deepvariant_R10 {
    tag "calling variants.."
    label 'deepvariant'
    memory '80 GB'
    time '24h'
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A30:1 --cpus-per-task=24 --qos=bonus'
    publishDir "$params.outdir/deepvariant", mode: 'copy'
    module 'singularity/3.7.4'
    module 'samtools/1.19.2'

    input:
        tuple val(lib), path(ontfile)
    output:
        tuple val(lib), path("*.vcf.gz"), emit: gz
        tuple val(lib), path("*.vcf.gz.tbi"), emit: tbi
        tuple val(lib), path("*.html"), emit: html
    script:
        def input_files = ("$ontfile".endsWith(".bam")) ? "${ontfile}" : ''
        """
        mkdir -p ${lib}_DeepVariant/intermediate_results_dir
	mkdir -p ~/vast_scratch/qg-rrms/cache
	mkdir -p ~/vast_scratch/qg-rrms/tmp

        export SINGULARITY_CACHEDIR=~/vast_scratch/qg-rrms/cache/
        export SINGULARITY_TMPDIR=~/vast_scratch/qg-rrms/tmp/
        export SINGULARITY_LOCALCACHEDIR=~/vast_scratch/qg-rrms/cache/

        BIN_VERSION="1.5.0"

        samtools index -@ 8 ${lib}_sup_5mCG_5hmCG.CHM13v2.bam

        singularity run --nv -B /vast -B /stornext -B /wehisan \
            docker://google/deepvariant:"\$BIN_VERSION-gpu" \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type=ONT_R104 \
            --ref="$params.fasta" \
            --reads=$input_files \
            --regions "$params.targets_bed" \
            --output_vcf="${lib}_output.vcf.gz" \
            --intermediate_results_dir "${lib}_DeepVariant/intermediate_results_dir" \
            --num_shards=24
        """


}

process filter_vcf {
    tag "filtering variants for PASS tag.."
    label 'filter'
    memory '80 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    clusterOptions '--qos=bonus'
    publishDir "$params.outdir/deepvariant", mode: 'copy'
    module 'bcftools/1.17'
    module 'htslib/1.17'
    module 'samtools/1.19.2'

    input:
        tuple val(lib), path(deepvariant_vcf)

    output:
        tuple val(lib), path("*_PASS.vcf.gz"), emit: gz
        tuple val(lib), path("*_PASS.vcf.gz.tbi"), emit: tbi

    script:
        """
        bcftools view -f PASS $deepvariant_vcf > "${lib}_PASS.vcf"
        bgzip "${lib}_PASS.vcf"
        tabix -p vcf "${lib}_PASS.vcf.gz"
        """
}

process phase_vcf {
    tag "phasing variants.."
    label 'phasing'
    memory '80 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    clusterOptions '--qos=bonus'
    publishDir "$params.outdir/deepvariant", mode: 'copy'
    module 'samtools/1.19.2'

    conda (params.enable_conda ? 'bioconda::whatshap=2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.1--py39h1f90b4d_0' :
        'quay.io/biocontainers/whatshap:2.1--py39h1f90b4d_0' }"

    input:
    tuple val(lib), path(deepvariant_vcf)
    tuple val(lib), path(deepvariant_idx)
    tuple val(lib), path(ontfile_bam)
    tuple val(lib), path(ontfile_idx)

    output:
    tuple val(lib), path("*.vcf.gz"), emit: gz
    tuple val(lib), path("*.vcf.gz.tbi"), emit: tbi

    script:
    """

    whatshap phase \
        --ignore-read-groups \
        --output ./${lib}.phased.vcf.gz \
        --reference "$params.fasta" \
        $deepvariant_vcf \
        $ontfile_bam

    tabix -p vcf ${lib}.phased.vcf.gz

    """
}

process phase_stats {
    tag "Making phasing annotation.."
    label 'annotation'
    memory '80 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    clusterOptions '--qos=bonus'
    publishDir "$params.outdir/deepvariant", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::whatshap=2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.1--py39h1f90b4d_0' :
        'quay.io/biocontainers/whatshap:2.1--py39h1f90b4d_0' }"

    input:
    tuple val(lib), path(deepvariant_vcf)
    tuple val(lib), path(deepvariant_idx)

    output:
    tuple val(lib), path("*.gtf"), emit: gtf

    script:
    """
    whatshap stats \\
    --gtf="${lib}.phased.gtf" \\
    "${lib}.phased.vcf.gz" \

    """


}

process mosdepth {
    tag "Calculating coverage.."
    label 'coverage'
    memory '80 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    clusterOptions '--qos=bonus'
    publishDir "$params.outdir/mosdepth", mode: 'copy'
    module 'mosdepth/0.3.1'

    input:
    tuple val(lib), path(ontfile_bam)

    output:
    tuple val(lib), path("*.mosdepth.global.dist.txt"), emit: txt

    script:
    """
    mosdepth -t 8 -x -n --by $params.chromsizes ${ontfile_bam}
    """

}



process index_hpbam {
    tag "Indexing bam files.."
    label 'index'
    memory '80 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    clusterOptions '--qos=bonus'
    publishDir "$params.outdir/whatshap", mode: 'copy'
    module 'samtools/1.19.2'

    input:
    tuple val(lib), path(ontfile_bam)

    output:
    tuple val(lib), path("*.bam.bai"), emit: bai

    script:
    """
    samtools index -@ 8 $ontfile_bam
    """
}



/*

process modbam2bed {
   tag "Making bedfiles.."
    label 'modbam2bed'
    memory '80 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    clusterOptions '--cpus-per-task=10 --qos=bonus'
    publishDir "$params.outdir/modbambed", mode: 'copy'
    module 'samtools/1.19.2'
    module 'bcftools/1.17'
    module 'htslib/1.17'

    input:
    tuple val(lib), path(ontfile_bam)


    output:
    tuple val(lib), path("*.bed.gz"), emit: gz

    script:
    """
    samtools index -@ 8 $ontfile_bam

    #mCG
    $projectDir/modbam2bed/modbam2bed -e -m 5mC --cpg -t 10 $params.fasta $ontfile_bam | bgzip -c > "${lib}_CHM13v2.mCG.bed.gz"

    #repeat for hmC
    $projectDir/modbam2bed/modbam2bed -e -m 5hmC --cpg -t 10 $params.fasta $ontfile_bam | bgzip -c > "${lib}_CHM13v2.hmcpg.bed.gz"

    #mCG haplotypes 1 and 2
    #HP1
    $projectDir/modbam2bed/modbam2bed -e -m 5mC --cpg -t 10 --haplotype 1 $params.fasta $ontfile_bam | bgzip -c > "${lib}_CHM13v2.hp1.mCG.bed.gz"

    #HP2
    $projectDir/modbam2bed/modbam2bed -e -m 5mC --cpg -t 10 --haplotype 2 $params.fasta $ontfile_bam | bgzip -c > "${lib}_CHM13v2.hp2.mCG.bed.gz"

    #hmC for haplotypes 1 and 2
    #HP1
    $projectDir/modbam2bed/modbam2bed -e -m 5hmC --cpg -t 10 --haplotype 1 $params.fasta $ontfile_bam | bgzip -c > "${lib}_CHM13v2.hp1.hmcpg.bed.gz"

    #HP2
    $projectDir/modbam2bed/modbam2bed -e -m 5hmC --cpg -t 10 --haplotype 2 $params.fasta $ontfile_bam | bgzip -c > "${lib}_CHM13v2.hp2.hmcpg.bed.gz"

    """
}
 */

/*
process create_bigwigs {
    tag "Making bigwigs.."
    label 'bigwig'
    memory '120 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    publishDir "$params.outdir/bigwigs", mode: 'copy'
    module 'samtools/1.19.2'
    module 'bcftools/1.17'
    module 'htslib/1.17'

    input:
    tuple val(lib), path(modbambed_bed)

    output:
    tuple val(lib), path("*.bw"), emit: bw

    script:
    """
    #mCG
    bgzip -dc "${lib}_CHM13v2.mCG.bed.gz" | cut -f 1,2,3,11 > "${lib}_bedGraph.bg"
    sort -k1,1 -k2,2n "${lib}_bedGraph.bg" > "${lib}_sorted.bedGraph"

    $projectDir/bedGraphToBigWig "${lib}_sorted.bedGraph" $params.chromsizes "${lib}_mCG.bw"

    #hmCG
    bgzip -dc "${lib}_CHM13v2.hmcpg.bed.gz" | cut -f 1,2,3,11 > "${lib}_bedGraph.bg"
    sort -k1,1 -k2,2n "${lib}_bedGraph.bg" > "${lib}_sorted.bedGraph"

    $projectDir/bedGraphToBigWig "${lib}_sorted.bedGraph" $params.chromsizes "${lib}_hmCpG.bw"

    #Make haplotyped bigWigs
    #hp1
    #mCG
    bgzip -dc "${lib}_CHM13v2.hp1.mCG.bed.gz" | cut -f 1,2,3,11 > "${lib}_hp1.bedGraph.bg"
    sort -k1,1 -k2,2n "${lib}_hp1.bedGraph.bg" > "${lib}_hp1.sorted.bedGraph"

    $projectDir/bedGraphToBigWig "${lib}_hp1.sorted.bedGraph" $params.chromsizes "${lib}_CHM13v2.hp1.mCG.bw"

    #hp2
    #mCG
    bgzip -dc "${lib}_CHM13v2.hp2.mCG.bed.gz" | cut -f 1,2,3,11 > "${lib}_hp2.bedGraph.bg"
    sort -k1,1 -k2,2n "${lib}_hp2.bedGraph.bg" > "${lib}_hp2.sorted.bedGraph"

    $projectDir/bedGraphToBigWig "${lib}_hp2.sorted.bedGraph" $params.chromsizes "${lib}_CHM13v2.hp2.mCG.bw"


    #hp1
    #hmCG
    bgzip -dc "${lib}_CHM13v2.hp1.hmcpg.bed.gz" | cut -f 1,2,3,11 > "${lib}_CHM13v2.hp1.hmcpg.bedGraph.bg"
    sort -k1,1 -k2,2n "${lib}_CHM13v2.hp1.hmcpg.bedGraph.bg" > "${lib}_CHM13v2.hp1.hmcpg.sorted.bedGraph"

    $projectDir/bedGraphToBigWig "${lib}_CHM13v2.hp1.hmcpg.sorted.bedGraph" $params.chromsizes "${lib}_hmCpG.hp1.bw"

    #hp2
    #hmCG
    bgzip -dc "${lib}_CHM13v2.hp2.hmcpg.bed.gz" | cut -f 1,2,3,11 > "${lib}_CHM13v2.hp2.hmcpg.bedGraph.bg"
    sort -k1,1 -k2,2n "${lib}_CHM13v2.hp2.hmcpg.bedGraph.bg" > "${lib}_CHM13v2.hp2.hmcpg.sorted.bedGraph"

    $projectDir/bedGraphToBigWig "${lib}_CHM13v2.hp2.hmcpg.sorted.bedGraph" $params.chromsizes "${lib}_hmCpG.hp2.bw"

    """
}
 */

//
// WORKFLOW: Run main nf-core/rrms analysis pipeline
//
workflow QG_RRMS {
    // parse input sample sheet
    ch_samples = INPUT_CHECK(ch_input)
        .groupTuple() // group samples by individual
        .map{individual, samples, bams -> tuple([id: individual, samples: samples], bams)} // merge individual and samples into a variable with id and sampels attribute.
    
    ch_merged_bam = SAMTOOLS_MERGE(ch_samples)
    ch_samples_bai = SAMTOOLS_INDEX_SAMPLES(ch_samples.transpose()).groupTuple()
    ch_merged_bai = SAMTOOLS_INDEX_MERGED(ch_merged_bam)
    ch_merged_bam.view()
    ch_vcf = DEEPVARIANT(
        params.deepvariant_region, 
        params.deepvariant_model, 
        ch_merged_bam, 
        ch_merged_bai, 
        ch_reference, 
        ch_reference_index
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    QG_RRMS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
