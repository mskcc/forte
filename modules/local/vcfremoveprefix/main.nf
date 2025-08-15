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
    tuple val(meta), path(vcf, name: 'input/*')

    output:
    tuple val(meta), path("*.{vcf}"), emit: vcf
    path "versions.yml"                  , emit: versions

    script:
    def modified_vcf = vcf.fileName.name
    """
    


    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """


}
