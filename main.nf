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

include { NANOPLOT } from '../modules/nf-core/nanoplot/main.nf'

//process batch {
//    label "cpu"

//}

process guppy_mod_basecall {
    label 'gpu'
    cpus 8
    memory '40GB'
    time '48h'
    queue 'gpuq'
    executor 'slurm'
    clusterOptions '--gres=gpu:A30:1 --cpus-per-task=8'
    publishDir "$params.outdir/sup_m5CG_basecalls", mode: 'copy'

    input:
        path fast5_dir
        path fasta
    output:
        tuple val(lib), path("*.fq.gz"), emit: fastq
        tuple val(lib), path("*.bam"), emit: bam
        tuple val(lib), path("*.txt.gz"), emit: summary
    script:
        """
        ~/Programs/ont-guppy_6.3.4/ont-guppy/bin/guppy_basecaller -i ${fast5_dir} -s sup_m5CG_basecalls --compress_fastq --config dna_r9.4.1_e8.1_modbases_5mc_cg_sup_prom.cfg --device "cuda:0" --recursive --min_qscore 7 --chunks_per_runner 768 --bam_out --align_ref ${fasta} --index
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

//process collect_batches {
//   label "cpu"

//}


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
        cat ${params.deepvariant_targets} | xargs -n1 -I% sh -c '/usr/local/bin/run_pepper_margin_deepvariant call_variant -b ${alignment_bam} -f ${params.reference} -o output/ --phased_output -p ${split_name}_%_PMDV -t 12 -g --ont_r9_guppy5_sup -r % || true'
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

    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    SAMPLESHEET_CHECK (ch_input)
    .csv
    .map {
        meta, fast5_dir ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fast5_dir }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)



    ch_guppy = guppy_mod_basecall(ch_fast5_dir,params.fasta)
    ch_nanoplot = NANOPLOT(ch_guppy.out.summary)
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
