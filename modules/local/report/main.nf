process REPORT_INDIVIDUAL {

    label "process_single"
    module "R/openBLAS/4.4.0"
    executor "local"

    input:
    tuple val(meta), path(mosdepth_html)

    output:
    path("${meta.id}_report.qmd")
    path(mosdepth_html)

    script:
    """
    echo '---
    title: "'"${meta.id}"'"
    format:
        html:
            grid:
                body-width: 1000px
            toc: true
    ---

    ## mosdepth coverage plots
    
    ```{=html}
    <figure style="width:100%;height:100%;">
        <iframe src="${mosdepth_html}" style="width:100%;height:80vh;overflow:hidden;margin:0px;padding:0px;border:none;"></iframe>
        <figcaption>Mosdepth coverage plots. "null" corresponds to the merged samples.</figcaption>
    </figure>
    ```
    ' > "${meta.id}"_report.qmd
    """

}

process REPORT_BOOK {

    label "process_single"
    module "R/openBLAS/4.4.0:quarto/1.3.450"
    executor "local"

    input:
    path(qmds)
    path(mosdepth_htmls)

    script:
    """
    # construct _quarto.yml
    echo 'project:
      type: book

    format:
      html:
        theme: cosmo
        number-depth: 3
    
    book: 
      title: "SkewX Report"
      author: "Quentin Gouil"
      date: "01/01/01"
      chapters:
        - index.qmd' > _quarto.yml
    for rep in *.qmd
    do
        echo "    - \$rep" >> _quarto.yml
    done

    # create an index.qmd
    echo "---
    title: Acknowledgements
    ---
    
    * My cat
    " >> index.qmd

    quarto render
    """


}