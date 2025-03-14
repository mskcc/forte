process FINGERPRINTPARSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.23.0--py39hdd5828d_0':
        'biocontainers/pysam:0.23.0--py39hdd5828d_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tsv"), emit: fp_txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_fingerprint_vcf.py \\
        --input ${vcf} \\
        --output ${prefix}.tsv \\
        --samplename ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fingerprintparse: 0.1.0 
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fingerprintparse: 0.1.0
    END_VERSIONS
    """
}
