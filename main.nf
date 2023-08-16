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
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input sample sheet not specified!' }

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

include { RRMS } from './workflows/rrms'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define processes and modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { NANOPLOT } from './modules/nf-core/nanoplot/main.nf'
include { NANOMETHVIZ } from './modules/nf-core/nanomethviz/main.nf'
include { SNIFFLES2 } from './modules/nf-core/sniffles/main.nf'
include { NANOCOMP } from './modules/nf-core/nanocomp/main.nf'
include { WHATSHAP } from './modules/nf-core/whatshap/main.nf'
include {INPUT_CHECK} from './subworkflows/local/input_check.nf'

//process batch {
//    label "cpu"

//}

process guppy_mod_basecall {
    debug true
    label 'gpu'
    cpus 8
    memory '40GB'
    time '48h'
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A30:1 --cpus-per-task=8'
    publishDir "$params.outdir/sup_m5CG_basecalls", mode: 'copy'

    input:
        tuple val(lib), path(raw_dir)
        //path fasta
    output:
        tuple val(lib), path("*.fq.gz"), emit: fastq
        tuple val(lib), path("*.bam"), emit: bam
        tuple val(lib), path("*.txt.gz"), emit: summary
    script:
        """
        ~/Programs/ont-guppy_6.3.4/ont-guppy/bin/guppy_basecaller -i ${raw_dir} -s sup_m5CG_basecalls --compress_fastq --config "$params.guppy_config" --device "cuda:0" --recursive --min_qscore 7 --chunks_per_runner 768 --bam_out --align_ref "$params.fasta" --index
        cd sup_m5CG_basecalls
        module load samtools/1.15

        samtools merge -@ 8 ${lib}.sup5mCG_passq7.bam pass/*.bam 
        samtools index ${lib}.sup5mCG_passq7.bam

        samtools merge -@ 8 ${lib}.sup5mCG_failq7.bam fail/*.bam 
        samtools index ${lib}.sup5mCG_failq7.bam

        cat pass/*.fastq.gz > ${lib}.sup.pass.fq.gz 
        cat fail/*.fastq.gz > ${lib}.sup.fail.fq.gz 

        mkdir -p logs
        mv *.log logs/
        tar czf ${lib}.log.tar.gz logs/ --remove-files
        gzip sequencing_summary.txt
        mv sequencing_summary.txt.gz ${lib}.sequencin_summary.txt.gz
        rm -r fail pass
        """
}

process dorado_mod_basecall {
    debug true
    label 'gpu'
    cpus 20
    memory '80GB'
    time '48h'
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A30:4 --cpus-per-task=20 --qos=bonus'
    publishDir "$params.outdir/sup_5mCG_5hmCG_alignments", mode: 'copy'
    module 'samtools/1.17'
    module 'dorado/0.3.2' 


    input:
        tuple val(lib), path(raw_dir)
    output:
        tuple val(lib), path("*.bam"), emit: bam
        tuple val(lib), path("*.bam.bai"), emit: bai
    script:
        """
        dorado basecaller "$params.dorado_config" ${raw_dir} --reference "$params.fasta" --modified-bases 5mCG_5hmCG | samtools sort > ${lib}_sup_5mCG_5hmCG.CHM13v2.bam
        samtools index ${lib}_sup_5mCG_5hmCG.CHM13v2.bam

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
    module 'samtools/1.18'

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

        export SINGULARITY_CACHEDIR=~/vast_scratch/qg-rrms/cache/
        export SINGULARITY_TMPDIR=~/vast_scratch/qg-rrms/tmp/
        export SINGULARITY_LOCALCACHEDIR=~/vast_scratch/qg-rrms/cache/

        BIN_VERSION="1.5.0"

        samtools index ${lib}_sup_5mCG_5hmCG.CHM13v2.bam

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
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A30:1 --cpus-per-task=24 --qos=bonus'
    publishDir "$params.outdir/deepvariant", mode: 'copy'
    module 'bcftools/1.17'
    module 'htslib/1.17'
    module 'samtools/1.18'

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
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A30:1 --cpus-per-task=24 --qos=bonus'
    publishDir "$params.outdir/deepvariant", mode: 'copy'

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
    tuple val(lib), path("*.vcf.gz"), emit: gz
    tuple val(lib), path("*.vcf.gz.tbi"), emit: tbi

    script:
    """
    whatshap phase \\
        --ignore-read-groups \\
        --output ./${lib}.phased.vcf.gz \\
        --reference "$params.fasta" \\
        $deepvariant_vcf \\
        $ontfile_bam \

    tabix -p vcf ${lib}.phased.vcf.gz

    """
}



//process collect_batches {
//   label "cpu"

//}

//TO DO: add R9.4.1 option to Pepper. --ont_r9_guppy5_sup

process run_pepper_margin_deepvariant {
    label 'deepvariant'
    cpus 4
    memory '64 GB'
    time '24h'
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A30:1 --nice'
    publishDir params.output_root, mode: 'copy', pattern: '*.vcf.gz', saveAs: {filename-> filename.substring(0, filename.indexOf("."))+"/"+filename}
    publishDir params.archive_root, mode: 'copy', pattern: '*PMDV*', saveAs: {filename-> filename.substring(0, filename.indexOf("."))+"/"+filename}

    input:
        path(input_fastq)
        val(sample_name)
        tuple path(alignment_bam), path(alignment_index)
    output:
        val(sample_name)
        path(alignment_bam)
        path(alignment_index)
        path("*PMDV.vcf.gz")
        path("*phased.vcf.gz"), emit: variant_calls
        path("*.tar.gz")
    script:
        split_name = sample_name+".deepvariant"
        """
        hostname
        set -euxo pipefail
        # Set up output directory
        mkdir output/
        # Run PEPPER-Margin-DeepVariant
        cat ${params.deepvariant_targets} | xargs -n1 -I% sh -c '/usr/local/bin/run_pepper_margin_deepvariant call_variant -b ${alignment_bam} -f ${fasta} -o output/ --phased_output -p ${split_name}_%_PMDV -t 12 -g --ont_r10_q20 -r % || true'
        find output/ -name "*phased*.vcf.gz" | xargs bcftools concat -o out.vcf
        bcftools sort -o ${sample_name}.PMDV_phased.vcf.gz -Oz out.vcf
        rm out.vcf
        find output/ -name "*.vcf.gz" | grep -v phase | xargs bcftools concat -a -o out.vcf
        bcftools sort -o ${sample_name}.PMDV.vcf.gz -Oz out.vcf
        cd output/
        tar czf ${sample_name}.PMDV_html_reports.tar.gz *.html
        tar czf ${sample_name}.PMDV_haplotagged_bams.tar.gz *.bam
        mv *.tar.gz ../
        cd ../
        """
}



process methylartist {
    label 'regular'
    cpus 8
    memory '32 GB'
    time '48h'
    queue 'regular'
    executor 'slurm'
    publishDir params.output_root, mode: 'copy', pattern: '*.zip', saveAs: {filename-> filename.substring(0, filename.indexOf("."))+"/"+filename}
    publishDir params.archive_root, mode: 'copy', pattern: '*.zip', saveAs: {filename-> filename.substring(0, filename.indexOf("."))+"/"+filename}

    input:
        val(sample_name)
        path(haplotag_bam)
        path(haplotag_bam_idx)
    output:
        file("*.zip")
        val(sample_name)
        tuple path(haplotag_bam), path(haplotag_bam_idx)
    shell:
        '''
        ln -s !{haplotag_bam} !{haplotag_bam}.bam
        ln -s !{haplotag_bam}.bai !{haplotag_bam}.bam.bai
        while read target_name target_plotregion target_highlightregion; do
            methylartist locus -b !{haplotag_bam}.bam -i ${target_plotregion} -l ${target_highlightregion} -g !{params.gff_file} --ref !{params.reference} --motif CG --labelgenes || true
            methylartist locus -b !{haplotag_bam}.bam -i ${target_plotregion} -l ${target_highlightregion} -g !{params.gff_file} --ref !{params.reference} --motif CG --labelgenes --phased || true
            methylartist region -p 8 -b !{haplotag_bam}.bam -i ${target_plotregion} -l ${target_highlightregion} -g !{params.gff_file} --ref !{params.reference} --motif CG --labelgenes || true
            methylartist region -p 8 -b !{haplotag_bam}.bam -i ${target_plotregion} -l ${target_highlightregion} -g !{params.gff_file} --ref !{params.reference} --motif CG --phased --labelgenes || true
        done < !{params.methylartist_regions}
        mkdir -p !{sample_name}.modplots/locus/phased
        mkdir -p !{sample_name}.modplots/locus/unphased
        mkdir -p !{sample_name}.modplots/region/phased
        mkdir -p !{sample_name}.modplots/region/unphased
        find . -maxdepth 1 -name "*phased*locus*" | xargs -n1 -I% mv % !{sample_name}.modplots/locus/phased
        find . -maxdepth 1 -name "*locus*" | xargs -n1 -I% mv % !{sample_name}.modplots/locus/unphased
        find . -maxdepth 1 -name "*phased*region*" | xargs -n1 -I% mv % !{sample_name}.modplots/region/phased
        find . -maxdepth 1 -name "*locus*" | xargs -n1 -I% mv % !{sample_name}.modplots/locus/unphased
        zip -r !{sample_name}.modplots.zip !{sample_name}.modplots/
        rm -rf !{sample_name}.modplots/
        unlink !{haplotag_bam}.bam
        unlink !{haplotag_bam}.bam.bai
        '''
}


process whatshap_haplotypes {
    label 'regular'
    cpus 12
    memory '24 GB'
    time '24h'
    queue 'regular'
    executor 'slurm'
    input:
        tuple file(fastqs), path(input_bam), path(input_bamidx), path(logpath), val(sample_name)
        path(phased_variant_calls)
    output:
        val(sample_name)
        file("*.bam")
        file("*.bai")
    script:
        """
        hostname
        tabix -p vcf ${phased_variant_calls}
        whatshap haplotag -o ${sample_name}.haplotagged.bam --reference ${params.reference} ${phased_variant_calls} ${input_bam} --ignore-read-groups
        samtools index ${sample_name}.haplotagged.modbam
        """
}




//
// WORKFLOW: Run main nf-core/rrms analysis pipeline
//
workflow QG_RRMS {
    ch_sample = INPUT_CHECK(ch_input)
    ch_sample.view()
    //ch_sample = Channel.fromPath(params.input).splitCsv( header:true, sep:',' )
    //    .map {  row -> [row[0], row[1]] } // lib, fast5_dir
    //ch_sample.view()
    //ch_guppy = guppy_mod_basecall(ch_sample)
    //ch_nanoplot = NANOPLOT(ch_dorado.summary)
    ch_dorado = dorado_mod_basecall(ch_sample)
    ch_nanocomp = NANOCOMP(ch_dorado.bam)
    ch_deepvariant = deepvariant_R10(ch_dorado.bam)
    ch_filter = filter_vcf(ch_deepvariant.gz)
    ch_phasing = phase_vcf(ch_filter.gz, ch_filter.tbi, ch_dorado.bam, ch_dorado.bai)
    ch_whatshap = WHATSHAP(ch_phasing.gz, ch_phasing.tbi, ch_dorado.bam, ch_dorado.bai)
    ch_sniffles = SNIFFLES2(ch_dorado.bam, ch_dorado.bai, params.tr_bed)
    ch_nanomethviz = NANOMETHVIZ(ch_dorado.bam, ch_dorado.bai)
    
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
