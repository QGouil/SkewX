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
        .map { create_fast5_dir_channel(it) }
        .set { ch_sample }

    emit:
    ch_sample                                     // channel: [ val(lib), [ fast5_dir ] ]
}

// Function to get list of [ lib, [ fast5_dir ] ]
def create_fast5_dir_channel(LinkedHashMap row) {
    // create lib map
    def lib = [:]
    lib.id         = row.sample
    lib.fast5_dir  = row.fast5_dir

    // add path(s) of the fastq file(s) to the lib map
    def fast5_lib = []
        fast5_lib = [ lib.id, [ path(lib.fast5_dir) ] ]
    return fast5_lib
}
