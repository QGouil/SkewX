process REPORT_INDIVIDUAL {

    tag "$meta.id"
    label "process_single"

    input:
    tuple val(meta), path(htmls), path(whatshap_stats_blocks), path(clustered_reads_tsv), path(skew_tsv), path(report_template)

    output:
    path("${meta.id}_report.qmd"), emit: qmds
    path(htmls), emit: htmls
    path(whatshap_stats_blocks), emit: whatshap_blocks
    path(clustered_reads_tsv), emit: clustered_reads

    script:
    """
    # copy individual template
    cp "${report_template}" "${meta.id}_report.qmd"

    # substitute individual id into report
    sed -i "s/meta_id/${meta.id}/g" "${meta.id}_report.qmd"

    # sub whatshap stats blocks file path into report
    sed -i "s/blocks_stats_file/${whatshap_stats_blocks}/g" "${meta.id}_report.qmd"
    """

}

process REPORT_BOOK {

    label "process_single"
    publishDir "${params.outdir}", mode: "copy"
    conda "${moduleDir}/../R/environment.yml"

    input:
    path(book_template_files)
    path(qmds)
    path(mosdepth_htmls) 
    path(whatshap_blocks)
    path(clustered_reads)
    path(cgi_bed) 

    output:
    path("_book")

    script:
    """
    # initialize _quarto.yml
    cp _quarto_template.yml _quarto.yml
    
    # append chapters (each patient) to _quarto.yml
    for rep in *_report.qmd
    do
        echo "    - \$rep" >> _quarto.yml
    done

    # substitute path to CGI bed file into each report
    sed -i "s/CGI_bed_file/${cgi_bed}/g" *_report.qmd

    quarto render
    """


}