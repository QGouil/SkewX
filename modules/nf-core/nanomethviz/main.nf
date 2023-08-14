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
                group = "5mCG"
            ),
            mod_code = "m"
    )
    p <- NanoMethViz:::plot_region(mbr, "chr15", 21626770, 21627383, heatmap=TRUE)
    svg("${lib}.MAGEL2.methplot.svg")
    print(p)
    dev.off()

    p1 <- NanoMethViz:::plot_region(mbr, "chr15", 21666002, 21666332, heatmap=FALSE)
    svg("${lib}.NDN.methplot.svg")
    print(p1)
    dev.off()

    p2 <- NanoMethViz:::plot_region(mbr, "chr15", 22936731, 22937323, heatmap=FALSE)
    svg("${lib}.SNRPN.methplot.svg")
    print(p2)
    dev.off()
    """
}
