process CFF_ANNOTATE {
    tag "$meta.id"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' :
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta),  path(cff), path(oncokb), path(agfusion)
    path(transcript_allowlist)

    output:
    tuple val(meta), path("${prefix}.unfiltered.cff"), emit: unfiltered_cff , optional: true
    tuple val(meta), path("${prefix}.final.cff")     , emit: filtered_cff , optional: true
    tuple val(meta), path("*.expanded_agfusion_transcripts.tsv"), emit: expanded_agfusion_transcripts, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${cff}"
    def oncokb_param = oncokb ? "--oncokb ${oncokb}" : ""
    """
    add_annotations_cff.R \\
        --cff ${cff} \\
        ${oncokb_param} \\
        --agfusion ${agfusion} \\
        --out-prefix ${prefix} \\
        --transcript_allowlist ${transcript_allowlist} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        add_annotations_cff.R: 0.0.2
    END_VERSIONS
    """
}
