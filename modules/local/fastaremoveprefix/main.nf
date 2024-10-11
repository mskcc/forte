process FASTAREMOVEPREFIX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    when:
    task.ext.when == null || task.ext.when

    input:
    tuple val(meta), path(fasta, name: 'input/*')

    output:
    tuple val(meta), path("*.{fa,fasta}"), emit: fasta
    path "versions.yml"                  , emit: versions

    script:
    def modified_fasta = fasta.fileName.name
    """
    cat ${fasta} | sed "s/^>chr/>/g" | sed "s/^>M />MT /g" > ${modified_fasta}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """


}
