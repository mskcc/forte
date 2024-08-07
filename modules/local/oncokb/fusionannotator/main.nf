process ONCOKB_FUSIONANNOTATOR {
    tag "$meta.id"
    label 'process_low'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    //conda "shahcompbio::oncokb-annotator=2.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/oncokb-annotator:0.0.1' :
        'docker://cmopipeline/oncokb-annotator:0.0.1' }"

    input:
    tuple val(meta), path(cff)

    output:
    tuple val(meta), path("*.oncokb.tsv")        , emit: oncokb_fusions
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'BEGIN { FS = "\\t"; OFS = "\\t"} {print \$31,\$19"-"\$21}' ${cff} | tail -n+2 > ${cff}.reformat
    echo -e "Tumor_Sample_Barcode\tFusion" | cat - ${cff}.reformat  > ${cff}.1reformat
    FusionAnnotator.py \\
        -i ${cff}.1reformat \\
        -o ${prefix}.oncokb.tsv \\
        ${args}

    GenerateReadMe.py \\
        -o README \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncokb_annotator: \$(grep "OncoKB data version" README | cut -f 2- -d":" | cut -f 1 -d, | tr -d " ")
        oncokb_api: \$(grep "OncoKB API URL" README | cut -f 2- -d":" | tr -d " ")
    END_VERSIONS
    """
}
