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
	    .splitCsv( header:true, sep:',' )
        .map{row -> tuple(row.individual, row.tissue, row.modbam_5mCG_dir)}
        .set{ch_sample}
 
    emit:
    ch_sample                                     // channel: [ val(lib), path(fast5_dir) ]
}

