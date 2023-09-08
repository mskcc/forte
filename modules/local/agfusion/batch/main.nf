process AGFUSION_BATCH {
    tag "$meta.id"
    label 'process_low'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda 'bioconda::agfusion=1.252'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/agfusion:0.0.6' :
        'cmopipeline/agfusion:0.0.6' }"

    input:
    tuple val(meta), path(fusions)
    path(agfusion_db)
    path(pyensembl_cache)

    output:
    tuple val(meta), path("${prefix}/**")                    , emit: fusions_annotated
    tuple val(meta), path("${prefix}.fusion_transcripts.csv"), emit: fusion_transcripts_csv
    tuple val(meta), path("${prefix}.fusion_transcripts.tsv"), emit: fusion_transcripts_tsv
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export PYENSEMBL_CACHE_DIR=\$PWD/${pyensembl_cache}

    awk -F"\\t" 'NR == 1 && \$1 ~ /gene5_chr/ {next;}{print}' ${fusions} > ${fusions}.no_header

    agfusion batch \\
        -f ${fusions}.no_header \\
        -db ${agfusion_db} \\
        -o ${prefix} \\
        ${args}

    awk -F"," 'NR != 1 && FNR == 1 {next;}{print}' ${prefix}/*/*.fusion_transcripts.csv > ${prefix}.fusion_transcripts.csv
    cat ${prefix}.fusion_transcripts.csv | tr "," "\\t" > ${prefix}.fusion_transcripts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agfusion: \$(agfusion -v) (fork)
    END_VERSIONS
    """
}
