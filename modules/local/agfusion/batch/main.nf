process AGFUSION_BATCH {
    tag "$meta.id"
    label 'process_low'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda 'bioconda::agfusion=1.252'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/agfusion:0.0.3' :
        'cmopipeline/agfusion:0.0.3' }"

    input:
    tuple val(meta), path(fusions)
    path(agfusion_db)
    path(pyensembl_cache)

    output:
    tuple val(meta), path("${prefix}")                       , emit: fusions_annotated
    tuple val(meta), path("${prefix}.fusion_transcripts.csv"), emit: fusion_transcripts
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export PYENSEMBL_CACHE_DIR=\$PWD/${pyensembl_cache}

    agfusion batch \\
        -f ${fusions} \\
        -db ${agfusion_db} \\
        -o ${prefix} \\
        ${args}

    cat ${prefix}/*/*.fusion_transcripts.csv | awk -F"," -v OFS="\\t" 'NR != 1 && FNR == 1 {next;}{print}' > fusion_transcripts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agfusion: \$(agfusion -v)
    END_VERSIONS
    """
}
