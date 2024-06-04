process REPORT_INDIVIDUAL {

    tag "$meta.id"
    label "process_single"

    input:
    tuple val(meta), 
          path(htmls),
          path(whatshap_stats),
          path(whatshap_blocks), 
          path(clustered_reads_tsv), 
          path(skew_tsv), 
          path(cgi_bed), 
          path(report_template)

    output:
    path("${meta.id}_report.qmd"), emit: qmds
    path(htmls), emit: htmls
    path("_${whatshap_stats.baseName}.qmd"), emit: whatshap_stats
    path(whatshap_blocks), emit: whatshap_blocks
    path(clustered_reads_tsv), emit: clustered_reads
    path(skew_tsv), emit: skew_tsv

    script:
    """
    # copy individual template
    cp "${report_template}" "${meta.id}_report.qmd"

    # substitute individual id into report
    sed -i "s/ext_meta_id/${meta.id}/g" "${meta.id}_report.qmd"

    # sub whatshap stats blocks file path into report
    sed -i "s/ext_blocks_stats_file/${whatshap_blocks}/g" "${meta.id}_report.qmd"

    # sub path to CGI bed file into each report
    sed -i "s/ext_CGI_bed_file/${cgi_bed}/g" "${meta.id}_report.qmd"

    # sub tissue names into report
    sed -i 's/ext_all_tissues_list/${meta.sample.findAll { !(it instanceof List)}.join('", "')}/g' "${meta.id}_report.qmd"

    # turn text files into qmd for code formatting
    echo '```' | cat - ${whatshap_stats} > "_${whatshap_stats.baseName}.qmd"
    echo '```' >> "_${whatshap_stats.baseName}.qmd"
    """

}

process REPORT_BOOK {

    label "process_low"
    publishDir "${params.outdir}", mode: "copy"
    conda "${moduleDir}/../R/environment.yml"

    input:
    path(book_template_files)
    path(qmds)
    path(mosdepth_htmls) 
    path(whatshap_stats)
    path(whatshap_blocks)
    path(clustered_reads)
    path(skews)
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

    # add downloadthis quarto extension
    quarto add --no-prompt "$projectDir/assets/report-templates/downloadthis"

    quarto render
    """


}