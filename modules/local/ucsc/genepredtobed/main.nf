process UCSC_GENEPREDTOBED {
    tag '${meta.id}'
    label 'process_low'

    conda "bioconda::ucsc-genepredtobed=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-genepredtobed:377--ha8a8165_5':
        'quay.io/biocontainers/ucsc-genepredtobed:377--ha8a8165_5' }"

    input:
    tuple val(meta), path(genepred)

    output:
    tuple val(meta), path("*.bed")     , emit: bed
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    genePredToBed \\
        $args \\
        $genepred  \\
        ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377'
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
