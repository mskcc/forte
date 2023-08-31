process COUNT_FEATURES {
    tag "$meta.id"
    label "process_low"

/// must be using singularity 3.7+
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.60.0--r43ha9d7317_0' :
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.60.0--r43ha9d7317_0' }"

    input:
    tuple val(meta), path(abundance)
    path(gtf)

    output:
    tuple val(meta), path("*.kallisto.customsummary.txt"), emit: kallisto_count_feature
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample = "${meta.sample}"
    """
    count_features.R \\
        --abundance ${abundance} \\
        --gtf ${gtf} \\
        --sample ${sample}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        count_features.R: 0.0.1
    END_VERSIONS
    """

    stub:
    def sample = meta.sample ?: ''
    """
    touch ${sample}.kallisto.customsummary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        add_flags_and_cluster_information.R: 0.0.1
    END_VERSIONS
    """
}
