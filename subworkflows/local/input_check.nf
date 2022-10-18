//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map {  it -> return tuple(it[0],it[1]) } // lib, fast5_dir
        .set { ch_sample }

    emit:
    ch_sample                                     // channel: [ val(lib), [ fast5_dir ] ]
}

