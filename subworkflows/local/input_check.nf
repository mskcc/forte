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
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.fq_num     = row.fastq_num.trim().toInteger()
    meta.sample     = row.sample.trim()
    meta.id         = meta.fq_num == 1 ? meta.sample : row.fastq_pair_id.trim()
    meta.fastq_pair_id = row.fastq_pair_id
    meta.single_end = row.single_end ? row.single_end.toBoolean() : false
    meta.bait       = row.bait ? row.bait.trim() : ""
    if (meta.bait != "") {
        try {
            params.baits[meta.bait]
        } catch(Exception e){
            throw e
            exit 1, "ERROR: Target files for baitset ${meta.bait} not available. Please check the samplesheet and add correct target files if applicable."
        }
    }

    meta.umi        = row.umi ? row.umi.trim() : ""
    meta.umi2       = row.umi2 ? row.umi2.trim() : ""
    meta.has_umi    = meta.umi == "" ? false : true
    if (meta.umi == "" && meta.umi2 != ""){
        exit 1, "ERROR: Please check input samplesheet -> there cannot be a umi2 pattern if there is no umi pattern:\n${meta.umi} and ${meta.umi2}"
    }
    try {
        meta.umi  = meta.umi.toInteger() * "N"
    } catch(Exception e) { }
    try {
        meta.umi2 = meta.umi2.toInteger() * "N"
    } catch(Exception e) { }
    meta.strand     = row.strand ? (row.strand.trim() == "" ? "no" : row.strand.trim()) : "no"
    if (! ["yes","no","reverse"].contains(meta.strand)){
        exit 1, "ERROR: Please check input samplesheet -> strand value is invalid!\n${row.strand ? row.strand : ""}"
    }

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}
