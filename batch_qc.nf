include {
    MULTIQC as MULTIQC_COLLECT
} from './modules/nf-core/multiqc/main'
include { COUNT_FEATURES } from './modules/local/count_features/main'
include { GUNZIP as GUNZIP_GTF } from './modules/nf-core/gunzip/main'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_forte_pipeline'
params.gtf = getGenomeAttribute('gtf')

ch_multiqc_config = Channel.fromPath("$projectDir/assets/analysis_multiqc_config.yml", checkIfExists: true)

workflow{
    BATCH_QC()
}

workflow BATCH_QC{

    multiqc_files = Channel.fromPath(params.input)
        .splitText()
    .map{line -> line.trim()}
    multiqc_files_dedup = multiqc_files
    .map{ sample_folder ->
        def matching_paths = []
        def base_dir = new File(sample_folder.trim())
        base_dir.eachDirRecurse { dir ->
            if (dir.toString().matches(".*/analysis/.*/dedup")){ matching_paths << file(dir)}
            if (dir.toString().matches(".*/analysis/.*/kallisto")){
                if (file( "$dir/*log.txt")) { matching_paths << file( "$dir/*log.txt") }
                if (file( "$dir/*.kallisto.customsummary.txt")) { matching_paths << file( "$dir/*.kallisto.customsummary.txt") }
            }
        }
        //def kallistofiles = files(sample_folder + "{/**,}/kallisto/*log.txt")
        //kallistofiles += files(sample_folder + "{/**,}/kallisto/*.kallisto.customsummary.txt")

        [matching_paths].flatten()

    }

    QC_DEDUP(multiqc_files_dedup)

    multiqc_files_dup = multiqc_files
        .map{ sample_folder ->
        def matching_paths = []
        def base_dir = new File(sample_folder)
        base_dir.eachDirRecurse { dir ->
            if (dir.toString().matches(".*/analysis/.*/dup")){ matching_paths << file(dir)}
            if (dir.toString().matches(".*/analysis/.*/star")){ matching_paths << file(dir)}
            if (dir.toString().matches(".*/analysis/.*/fastp")){
                matching_paths << files( "$dir/*.fastp.json")
            }

        }

        [matching_paths].flatten()

        }

    QC_DUP(multiqc_files_dup)
}

workflow QC_DEDUP {

    take:
    multiqc_files_dedup


    main:

    MULTIQC_COLLECT(
        multiqc_files_dedup.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        [],
        [],
        [],
        []
    )
}


workflow QC_DUP {

    take:
    multiqc_files_dup

    main:
    MULTIQC_COLLECT(
        multiqc_files_dup.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        [],
        [],
        [],
        []
    )



}

