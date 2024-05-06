process REPORT_INDIVIDUAL {

    label "process_single"
    module "R/openBLAS/4.4.0"
    executor "local"

    input:
    tuple val(meta), path(htmls)

    output:
    path("${meta.id}_report.qmd")
    path(htmls)

    script:
    """
    echo '---
    title: "'"${meta.id}"'"
    format:
        html:
            grid:
                body-width: 900px
            toc: true
    ---

    ## QC

    ### mosdepth coverage plots
    
    ```{=html}
    <figure style="width:100%;height:100%;">
        <iframe src="${meta.id}.dist.html" style="width:900px;height:1000px;overflow:hidden;margin:0px;padding:0px;border:none;"></iframe>
        <figcaption>Mosdepth coverage plots.</figcaption>
    </figure>
    ```

    ### nanocomp plots

    ```{=html}
    <iframe src="${meta.id}_NanoComp_log_length_violin.html" style="width:100%;height:100vh;overflow:hidden;margin:0px;padding:0px;border:none;"></iframe>
    ```

    ```{=html}
    <iframe src="${meta.id}_NanoComp_number_of_reads.html" style="width:100%;height:100vh;overflow:hidden;margin:0px;padding:0px;border:none;"></iframe>
    ```

    ```{=html}
    <iframe src="${meta.id}_NanoComp_N50.html" style="width:100%;height:100vh;overflow:hidden;margin:0px;padding:0px;border:none;"></iframe>
    ```

    ```{=html}
    <iframe src="${meta.id}_NanoComp_percentIdentity_violin.html" style="width:100%;height:100vh;overflow:hidden;margin:0px;padding:0px;border:none;"></iframe>
    ```

    ```{=html}
    <iframe src="${meta.id}_NanoComp_quals_violin.html" style="width:100%;height:100vh;overflow:hidden;margin:0px;padding:0px;border:none;"></iframe>
    ```
    ' > "${meta.id}"_report.qmd
    """

}

process REPORT_BOOK {

    label "process_single"
    module "R/openBLAS/4.4.0:quarto/1.3.450"
    executor "local"
    publishDir "${params.outdir}", mode: "copy"

    input:
    path(qmds)
    path(mosdepth_htmls)

    output:
    path("_book")

    script:
    """
    # initialize _quarto.yml
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
    
    # append chapters (each patient) to _quarto.yml
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