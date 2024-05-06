process REPORT_INDIVIDUAL {

    label "process_single"
    module "R/openBLAS/4.4.0"
    executor "local"

    input:
    tuple val(meta), path(htmls), path(report_template)

    output:
    path("${meta.id}_report.qmd")
    path(htmls)

    script:
    """
    # copy individual template
    cp "${report_template}" "${meta.id}_report.qmd"

    # substitute individual id into report
    sed -i "s/meta_id/${meta.id}/g" "${meta.id}_report.qmd"
    """

}

process REPORT_BOOK {

    label "process_single"
    module "R/openBLAS/4.4.0:quarto/1.3.450"
    executor "local"
    publishDir "${params.outdir}", mode: "copy"

    input:
    path(book_template_files)
    path(qmds)
    path(mosdepth_htmls)    

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

    quarto render
    """


}