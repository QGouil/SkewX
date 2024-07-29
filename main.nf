#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SkewX: A Nextflow pipeline for skewed X inactivation analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/QGouil/SkewX

    Publication: https://doi.org/10.1101/2024.03.20.585856

----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


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
        it -> tuple(id: it.baseName, it, "${it}.fai")
    }
} else { exit 1, "Reference FASTA not specified!" }
if (params.cgi_bedfile) {
    ch_cgibed = Channel.fromPath(params.cgi_bedfile, checkIfExists: true)
        .map{
            it -> tuple(id: it.baseName, it)
        }
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
include {MINIMAP2} from "./modules/local/minimap2/main.nf"
include {SAMTOOLS_MERGE} from "./modules/local/samtools/merge/main.nf"
// two processes needed as two sets of bams are indexed.
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_SAMPLES} from "./modules/local/samtools/index/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED} from "./modules/local/samtools/index/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_HAPLOTAG} from "./modules/local/samtools/index/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_HAPLOTAG_MERGED} from "./modules/local/samtools/index/main.nf"
include {FILTER_PASS} from "./modules/local/bcftools/view_pass/main.nf"
include {WHATSHAP_PHASE} from "./modules/local/whatshap/phase/main.nf"
include {WHATSHAP_STATS} from "./modules/local/whatshap/stats/main.nf"
include {WHATSHAP_HAPLOTAG} from "./modules/local/whatshap/haplotag/main.nf"
include {WHATSHAP_HAPLOTAG as WHATSHAP_HAPLOTAG_MERGED} from "./modules/local/whatshap/haplotag/main.nf"
include {MOSDEPTH} from "./modules/local/mosdepth/main.nf"
include {MOSDEPTH as MOSDEPTH_MERGED} from "./modules/local/mosdepth/main.nf"
include {SAMTOOLS_VIEWHP} from "./modules/local/samtools/view_hp/main.nf"
include {R_CLUSTERBYMETH} from "./modules/local/R/cluster_by_meth/main.nf"
include {reporting} from "./subworkflows/reporting.nf"
include {separated_deepvariant} from "./subworkflows/local/deepvariant/main.nf"

//
// WORKFLOW: Run main SkewX analysis pipeline
//
workflow SKEWX {
    // parse input sample sheet
    ch_checked_input = INPUT_CHECK(ch_input)

    // put individual id and sample into first element of tuple.
    ch_separate_samples = ch_checked_input
        .map{individual, sample, bam -> tuple([id: individual, sample: sample], bam)}

    // if reads are not mapped, align with minimap2, otherwise assume input bams are aligned
    if (params.ubam) {
        ch_aligned = MINIMAP2(ch_separate_samples, ch_reference)
    } else {
        ch_aligned = ch_separate_samples
    }

    // index aligned bams
    ch_samples = SAMTOOLS_INDEX_SAMPLES(ch_aligned)

    // prepare for merging by grouping
    ch_grouped_bams = ch_samples
        .map{meta, bam, bai -> tuple(meta.id, meta.sample, bam, bai)} // unpack id, sample to enable grouping by individual id
        .groupTuple() // group samples by individual
        .map{individual, samples, bams, bais -> tuple([id: individual, sample: samples], bams, bais)} // merge individual and samples into a variable with id and samples attribute.

    // separate individuals with only one sample to bypass merging
    (ch_multiple_bams, ch_single_bams) = ch_grouped_bams
        .branch{
            multi: it[1].size()>1
            single: it[1].size()==1
        }
    ch_single_bams = ch_single_bams.map{
        meta, bams, bais -> tuple(meta, bams[0], bais[0])
    } // unpack bam from list container

    SAMTOOLS_MERGE(ch_multiple_bams)
        | SAMTOOLS_INDEX_MERGED
        | mix(ch_single_bams) // add back individuals with single bams
        | set {ch_merged_bam}

    // variant call merged bam with deepvariant
    // first, duplicate ch_reference for each merged bam
    ch_reference_rep_merged = ch_merged_bam
        .combine(ch_reference.collect())
        .map{meta, merged_bams, merged_bams_idx, meta_ref, ref, ref_idx -> tuple(meta_ref, ref, ref_idx)}
    dv_args = channel.from([
        regions: params.deepvariant_region,
        model_type: params.deepvariant_model,
        num_shards: params.deepvariant_num_shards
    ])
    (ch_vcf, ch_deepvariant_report) = separated_deepvariant(dv_args, ch_merged_bam, ch_reference)
        .multiMap{
            vcf: tuple(it[0], it[1], it[2], it[5])
            dv_report: tuple(it[0], it[7])
        }

    // filter variants by PASS
    ch_vcf_pass = FILTER_PASS(ch_vcf)

    // phase variants
    ch_vcf_phased = WHATSHAP_PHASE(
        ch_vcf_pass,
        ch_reference_rep_merged
    )

    ch_whatshap_stats_blocks = WHATSHAP_STATS(ch_vcf_phased.map{meta, bam, bam_idx, vcf, vcf_idx -> tuple(meta, vcf, vcf_idx)})

    // repeat phased vcf and reference channels, so there are enough items to
    // match the number of samples being haplotagged.
    (ch_tmp_samples, ch_reference_rep) = ch_samples
        .map{meta, bam, bam_idx -> tuple(meta.id, meta.sample, bam, bam_idx)} // unpack meta so vcfs can be added to each sample, based on id
        .combine( // combine vcfs based on individual id
            ch_vcf_phased.map{meta, merged_bam, merged_bam_idx, vcf, vcf_idx -> tuple(meta.id, meta.sample, vcf, vcf_idx)},
            by: 0
        )
        .combine(ch_reference.collect()) // repeat reference for each sample
        .multiMap{id, single_sample, bam, bai, samples, vcf, vcf_idx, ref_id, ref, ref_idx -> // unpack the now horrendously long item into separate channels
            samples: tuple([id: id, sample: single_sample], bam, bai, vcf, vcf_idx)
            ref: tuple(ref_idx, ref, ref_idx)
        }

    // haplotag individual sample bams and index each
    WHATSHAP_HAPLOTAG(ch_tmp_samples, ch_reference_rep)
        | SAMTOOLS_INDEX_HAPLOTAG
        | set {ch_samples_haplotag}
    (ch_mosdepth, ch_mosdepth_report_results) = MOSDEPTH(ch_samples_haplotag)

    // haplotag merged sample bams and index each
    WHATSHAP_HAPLOTAG_MERGED(ch_vcf_phased, ch_reference_rep)
        | SAMTOOLS_INDEX_HAPLOTAG_MERGED
        | set {ch_merged_samples_haplotag}

    (_, ch_mosdepth_report_results_merged) = MOSDEPTH_MERGED(ch_merged_samples_haplotag)

    // repeat cigx bed to match each haplotagged sample
    (ch_tmp_samples_haplotag, ch_cgibed_rep) = ch_samples_haplotag
        .combine(ch_cgibed.collect())
        .multiMap{it ->
            samples_haplotag: tuple(it[0], it[1], it[2])
            cgibed: tuple(it[3], it[4])
        }

    ch_hpreads = SAMTOOLS_VIEWHP(ch_tmp_samples_haplotag, ch_cgibed_rep)

    ch_clustered_reads = R_CLUSTERBYMETH(ch_hpreads, ch_cgibed_rep)

    // prepare inputs for plotting mosdepth results
    ch_mosdepth_report_results_merged
        .filter{it[0].sample.size()>1} // remove merged samples with only one constituent samplebi
        .map{it -> tuple(it[0].id, it[0].sample, it[1])} // extract individual id, sample, and *.global.dist.txt from channel
        .mix(ch_mosdepth_report_results.map{it -> tuple(it[0].id, it[0].sample, it[1])}) // mix in mosdepth results for single samples
        .groupTuple() // group by individual id
        .map{ it -> tuple([id: it[0], sample: it[1]], it[2])} // merge id and samples into tuple
        .set{ch_mosdepth_all_report_results}

    // nanocomp
    ch_merged_samples_haplotag
        .filter{it[0].sample.size()>1} // only include merged bams constisting of multiple samples
        .map{it -> tuple(it[0].id, it[0].sample, it[1], it[2])}
        .mix(ch_samples_haplotag.map{it->tuple(it[0].id, it[0].sample, it[1], it[2])})
        .groupTuple()
        .map{it -> tuple([id: it[0], sample: it[1]], it[2], it[3])}
        .set{ch_all_samples_haplotag}

    // put together report
    book = reporting(
        ch_mosdepth_all_report_results,
        ch_all_samples_haplotag,
        ch_whatshap_stats_blocks,
        ch_clustered_reads,
        ch_cgibed,
        "${projectDir}/assets/report-templates/_extensions"
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    SKEWX ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
