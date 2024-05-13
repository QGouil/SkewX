include {MOSDEPTH_PLOTDIST} from "../modules/local/mosdepth/plotdist/main.nf"
include {NANOCOMP} from "../modules/local/nanocomp/main.nf"
include {REPORT_INDIVIDUAL} from "../modules/local/report/main.nf"
include {REPORT_BOOK} from "../modules/local/report/main.nf"
include {CUMULATIVE_CGI_COVERAGE} from "../modules/local/R/plots/main.nf"

workflow reporting {
    take:
        mosdepth_report_results // channel containing by-sample and merged mosdepth report results
        haplotagged_samples // channel containing by-sample and merged sample haplotagged bams
        whatshap_stats_blocks // channel containing phased vcfs block stats
        clustered_reads // channel containining clustered reads and skew information 
        cgi_bed // single-item channel containing CGI bed file
    main:
        // prepare mosdepth coverage report
        // add script into channels
        ch_plot_dist_script = Channel.fromPath("https://raw.githubusercontent.com/brentp/mosdepth/v0.3.6/scripts/plot-dist.py")
        (ch_global_dist_bysample, ch_plot_dist_script_rep) = mosdepth_report_results
            .combine(ch_plot_dist_script)
            .multiMap{it ->
                dists: tuple(it[0], it[1])
                scripts: it[2]
            }
        ch_mosdepth_dist_report = MOSDEPTH_PLOTDIST(ch_plot_dist_script_rep, ch_global_dist_bysample)

        // prepare nanocomp reports
        ch_nanocomp = NANOCOMP(haplotagged_samples)

        // combine nanocomp + mosdepth reports
        ch_combined_qc_reports = ch_nanocomp
            .map{it -> tuple(it[0].id, it[0].sample, it[1])}
            .join(ch_mosdepth_dist_report.map{it -> tuple(it[0].id, it[0].sample, it[1])})
            .join(whatshap_stats_blocks.map{it -> tuple(it[0].id, it[0].sample, it[1], it[2])}) // blocks in 3rd element
            .join(clustered_reads.map{it -> tuple(it[0].id, it[0].sample, it[1], it[2])}.groupTuple())
            .map{it -> tuple([id: it[0], sample: it[1]], it[2] + [it[4]], it[6], it[7], it[9], it[10])}
            .combine(cgi_bed.map{it -> it[1]})
            .combine(channel.fromPath("${projectDir}/assets/report-templates/individual_report.qmd", checkIfExists: true))

        ch_reporting_files = REPORT_INDIVIDUAL(ch_combined_qc_reports)

        // create channel for templates
        ch_book_template_files = channel.fromPath([
            "${projectDir}/assets/report-templates/_quarto_template.yml",
            "${projectDir}/assets/report-templates/index.qmd"
        ], checkIfExists: true).collect() // collect to make sure all the files are one item
        book = REPORT_BOOK(
            ch_book_template_files, 
            ch_reporting_files.qmds.collect(), 
            ch_reporting_files.htmls.collect(),
            ch_reporting_files.whatshap_stats.collect(),
            ch_reporting_files.whatshap_blocks.collect(),
            ch_reporting_files.clustered_reads.collect(),
            cgi_bed.map{it -> it[1]}
        )
    emit:
        book
}