process ARRIBA_DOWNLOAD {
    tag { "${prefix}" }
    label 'process_low'

    conda "bioconda::arriba=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.3.0--haa8aa89_0' :
        'quay.io/biocontainers/arriba:2.3.0--haa8aa89_0' }"

    output:
    path "blacklist*"       , emit: blacklist
    path "protein_domains*" , emit: protein_domains
    path "known_fusions*"   , emit: known_fusions
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "hg37"

    """
    wget https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz -O arriba_v2.3.0.tar.gz
    tar -xzvf arriba_v2.3.0.tar.gz
    rm arriba_v2.3.0.tar.gz
    cp arriba_v2.3.0/database/blacklist_${prefix}_* .
    cp arriba_v2.3.0/database/known_fusions_${prefix}_* .
    cp arriba_v2.3.0/database/protein_domains_${prefix}_* .
    rm -r arriba_v2.3.0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_hg38_GRCh38_v2.3.0.tsv.gz
    touch protein_domains_hg38_GRCh38_v2.3.0.gff3
    touch cytobands_hg38_GRCh38_v2.3.0.tsv
    touch known_fusions_hg38_GRCh38_v2.3.0.tsv.gz
    touch protein_domains_hg38_GRCh38_v2.3.0.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
