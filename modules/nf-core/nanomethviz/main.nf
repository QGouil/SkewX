process NANOMETHVIZ {
    tag "Plotting regions of interest..."
    label "NanoMethViz"
    cpus 8
    memory '32 GB'
    time '48h'
    queue 'regular'
    executor 'slurm'
    publishDir "$params.outdir/NanoMethViz", mode: 'copy'
    module 'R/4.2.0'

    input:
    tuple val(lib), path(ontfile_bam)
    tuple val(lib), path(ontfile_idx)

    output:
        path "*.methplot.svg", emit: svg

    script:
    """
    #!/usr/bin/env Rscript

    #remotes::install_github('shians/nanomethviz')

    library(NanoMethViz)

    mbr <- ModBamResult(
            methy = ModBamFiles(
                samples = "${lib}",
                paths = "${ontfile_bam}"),
            samples = data.frame(
                sample = "${lib}",
                group = "group"
            ),
            mod_code = "m"
    )
    p <- NanoMethViz:::plot_region(mbr, "chr1", 247365298, 247375298, heatmap=TRUE)
    svg("test.methplot.svg")
    print(p)
    dev.off()
    """
}