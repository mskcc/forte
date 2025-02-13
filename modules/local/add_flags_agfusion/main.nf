process ADD_FLAGS_AGFUSION {
    tag "$meta.id"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' :
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta), path(cff), path(agfusion)
    path(transcript_allowlist)

    output:
    tuple val(meta), path("*.expanded_agfusion_transcripts.tsv"), emit: expanded_agfusion_transcripts
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.sample}"

    """
    add_flags_agfusion_clinical.R \\
        --cff $cff \\
        --agfusion $agfusion \\
        --transcript_allowlist $transcript_allowlist \\
        --out-prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        add_flags_agfusion_clinical.R: 0.0.1
    END_VERSIONS
    """
}
