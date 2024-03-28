process render_report {
    tag "creating report"
    label 'process_high'
    publishDir "${params.out_dir}/report/", mode: 'copy', pattern: "*.zip"
   // stageInMode 'copy'
    container "library://kzeglinski/nanologix/nanologix-report:v0.0.3"

    input:
    path(processed_tsv)
    path(sample_sheet)
    path(multiqc_plots)
    path(percentage_passing_trim_merge)
    path(template_dir)
    path(extensions_dir)
    path(quarto_base_yaml)
    path(qmd_templates)
    val(analysis_name)

    output:
    path('*.zip'), emit: report

    script:
    """
    #!/usr/bin/env bash
    tar -xvf _extensions.tar
    tar -xvf _template.tar

    quarto render --log quarto.log

    mkdir report_${analysis_name}
    cp -r _book/ report_${analysis_name}
    cp *.fasta report_${analysis_name}
    cp *.qmd report_${analysis_name}
    cp *.yml report_${analysis_name}
    cp *.log report_${analysis_name}
    cp mqc_fastqc_per_base_sequence_quality_plot_1.png report_${analysis_name}
    cp mqc_fastqc_per_sequence_quality_scores_plot_1.png report_${analysis_name}
    cp mqc_fastqc_per_sequence_gc_content_plot_Percentages.png report_${analysis_name}
    cp -r template/ report_${analysis_name}
    cp -r _extensions/ report_${analysis_name}
    zip -r report_${analysis_name}.zip report_${analysis_name}/
    """
}

process prepare_report_templates {
    tag "preparing report templates"
    label 'process_low'
    stageInMode 'copy'
    container "library://kzeglinski/nanologix/nanologix-report:v0.3.0"

    input:
    path(sample_sheet)
    path(quarto_base_yaml)
    path(qmd_templates)
    val(adapter_r1)
    val(adapter_r2)
    val(analysis_name)
    val(sequence_trim_5p)
    val(sequence_trim_3p)

    output:
    path('*.qmd'), emit: report_templates
    path('*.yml'), emit: edited_quarto_yaml

    script:
    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(vroom)

    # read in the sample sheet
    metadata <- vroom(fs::dir_ls(glob = "*.csv"),
        col_select = c("sample_num", "library", "antigen", "round", "replicate"))

    # create replicate ID if required
    if(!all(is.na(metadata[["replicate"]]))){
        metadata <- metadata %>%
            filter(round != 0) %>% # temporarily remove R0 samples
            group_by(library, antigen, round) %>%
            mutate(replicate_id = cur_group_id()) %>%
            # and give our replicates an informative name, not just a number
            mutate(replicate_id_informative =
                paste0(library, "_", antigen, "_round_", round)) %>%
            ungroup() %>%
            bind_rows(filter(metadata, round == 0))
    }

    # create panning ID if required
    if(!all(metadata[["round"]] == 0)){
        metadata <- metadata %>%
            filter(round != 0) %>% # temporarily remove R0 samples
            group_by(library, antigen) %>%
            mutate(panning_id = as.character(cur_group_id())) %>%
            ungroup() %>%
            bind_rows(filter(metadata, round == 0))

        metadata <- metadata %>%
            group_by(library) %>%
            summarise(panning_id = paste0(unique(na.omit(panning_id)), collapse = "_")) %>%
            mutate(round = 0, replicate = NA) %>%
            right_join(metadata, by = c("library", "round", "replicate")) %>%
            mutate(panning_id = coalesce(panning_id.x, panning_id.y)) %>%
            select(-panning_id.x, -panning_id.y)
    }

    # ones we're not changing
    readLines("template_index.qmd") %>%
        writeLines(con = "index.qmd")

    readLines("template_references.qmd") %>%
        writeLines(con = "references.qmd")

    readLines("template_sequencing_qc.qmd") %>%
        writeLines(con = "sequencing_qc.qmd")

    readLines("template_library_qc.qmd") %>%
        writeLines(con = "library_qc.qmd")

    readLines("template_panning.qmd") %>%
        writeLines(con = "panning.qmd")

    # introduction
    readLines("template_intro.qmd") %>%
        stringr::str_replace(pattern = "param_analysis_name", replace = "${analysis_name}") %>%
        stringr::str_replace(pattern = "param_adapter_r1", replace = "${adapter_r1}") %>%
        stringr::str_replace(pattern = "param_adapter_r2", replace = "${adapter_r2}") %>%
        writeLines(con = "intro.qmd")

    # the YAML file
    # need to add all of the individual pan files
    yaml_file <- yaml::read_yaml("template_quarto.yml")

    if(!all(metadata[["round"]] == 0)){
        all_pans <- metadata %>%
            filter(round != 0) %>%
            pull(panning_id) %>%
            unique() %>%
            paste0("pan_", ., ".qmd")
        yaml_file[["book"]][["chapters"]][[5]][["chapters"]] <- all_pans
    } else {
        # if we don't have any panning
        yaml_file[["book"]][["chapters"]][[5]] <- NULL
    }

    # need this to make sure it keeps the logicals as true and false not yes and no (quarto doesn't like that)

    yaml::write_yaml(yaml_file, "_quarto.yml" ,handlers = list(
    logical = function(x) {
      result <- ifelse(x, "true", "false")
      class(result) <- "verbatim"
      return(result)
    }))

    # now actually made these .qmd files
    pan_id <- all_pans %>% str_remove("pan_") %>% str_remove(".qmd")
    for(i in seq_along(all_pans)) {
        this_pan_id <- pan_id[i]
        this_pan_name <- metadata %>%
            filter(panning_id == this_pan_id) %>%
            slice(1) %>%
            summarise(pan_name = paste0(library, " vs ", antigen)) %>%
            pull()
        readLines("single_pan_template.qmd") %>%
            stringr::str_replace(pattern = "param_id_num", replace = this_pan_id) %>%
            stringr::str_replace(pattern = "param_name", replace = this_pan_name) %>%
            stringr::str_replace(pattern = "param_trim_5p", replace = "'${sequence_trim_5p}'") %>%
            stringr::str_replace(pattern = "param_trim_3p", replace = "'${sequence_trim_3p}'") %>%
            writeLines(all_pans[i])
    }
    """
}
