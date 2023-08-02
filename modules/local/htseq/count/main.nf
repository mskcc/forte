process HTSEQ_COUNT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::htseq=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py39h919a90d_0' :
        'quay.io/biocontainers/htseq:2.0.2--py39h919a90d_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path gff

    output:
    tuple val(meta), path("*.htseq.count.txt"), emit: counts
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    htseq-count \\
        $args \\
        $bam \\
        $gff > \\
        ${prefix}.htseq.count.txt

    sed -i '1{h;s/.*/'"${prefix}"'/;G}' ${prefix}.htseq.count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(htseq-count --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.htseq.count.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(htseq-count --version )
    END_VERSIONS
    """
}
