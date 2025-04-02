process FINGERPRINT_CORRELATION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-data.table:1.11.6--r341hc070d10_0':
        'biocontainers/r-data.table:1.11.6--r341hc070d10_0' }"

    input:
    path(fingerprints, name: 'input/*')

    output:
    path("correlation.tsv"), emit: correlation
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo -e "sample\\tpath" > input.tsv
    for i in input/* ; do
        echo -e "\$(basename \$i | cut -f 1 -d.)\\t\$i" >> input.tsv
    done

    fingerprint_correlation.R \\
        --input_table input.tsv \\
        --output correlation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fingerprint_correlation.R: 0.1.0
        r-data.table: 1.11.6
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch correlation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fingerprint_correlation.R: 0.1.0
        r-data.table: 1.11.6
    END_VERSIONS
    """
}
