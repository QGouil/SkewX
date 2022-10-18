//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv //contains the samplesheet_valid.csv
        .splitCsv( header:true, sep:',' ).view{ row -> "${row.sample} - ${row.fast5_dir}"}
        //.map{row -> tuple(row.sample, row.fast5_dir)} // lib, fast5_dir
        //.set{ch_sample}

    emit:
    ch_sample                                     // channel: [ val(lib), path(fast5_dir) ]
}

