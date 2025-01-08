//
// Check input samplesheet and get read channels
//

include { MAF_SAMPLESHEET_CHECK } from '../../modules/local/maf_samplesheet_check/maf_samplesheet_check'


workflow MAF_INPUT_CHECK {
    take:
    maf_samplesheet // file: /path/to/samplesheet.csv
    all_samples

    main:

    input = maf_samplesheet ? Channel.fromPath(maf_samplesheet) : Channel.empty()
    MAF_SAMPLESHEET_CHECK ( input )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_maf_channel(it) }
        .set { mafs }

    mafs.map{meta, maf -> [meta.sample,meta.sample]}
        .join(all_samples.map{it -> [it, it]}, remainder: true)
        .map{ it -> [it[0], it[2]] }
        .map{ maf_sample, sample->
            if (sample == null) {
                println "WARNING: Sample in the maf input sheet does not match samples in the regular input sheet:\n${maf_sample}"
            }
        }

    emit:
    mafs
    versions = MAF_SAMPLESHEET_CHECK.out.versions
}

def create_maf_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.sample = row.sample.trim()
    def maf_file = file(row.maf)
    if (!file(maf_file).exists()){
        exit 1, "ERROR: Please check fillout input samplesheet -> MAF file does not exist!\n${row.maf}"
    }
    return [meta,maf_file]
}
