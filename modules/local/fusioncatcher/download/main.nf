process FUSIONCATCHER_DOWNLOAD {
    tag 'fusioncatcher_download'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/clinicalgenomics/fusioncatcher:1.33' :
        'docker.io/clinicalgenomics/fusioncatcher:1.33' }"

    output:
    path "fusioncatcher_*", emit: reference
    path "versions.yml"   , emit: versions

    script:

    def args             = task.ext.args ?: ''
    def args2            = task.ext.args2 ?: ''
    def species          = 'homo_sapiens'
    def download_version = 'human_v102'
    def url = "http://sourceforge.net/projects/fusioncatcher/files/data/${download_version}.tar.gz.aa"
    """
    if wget --spider "$url" 2>/dev/null; then
        wget $args $url
        wget $args http://sourceforge.net/projects/fusioncatcher/files/data/${download_version}.tar.gz.ab
        wget $args http://sourceforge.net/projects/fusioncatcher/files/data/${download_version}.tar.gz.ac
        wget $args http://sourceforge.net/projects/fusioncatcher/files/data/${download_version}.tar.gz.ad
        cat ${download_version}.tar.gz.* | tar xz
        mv ${download_version} fusioncatcher_${download_version}
        rm ${download_version}.tar*
    else
        fusioncatcher-build.py \\
            -g ${species} \\
            -o fusioncatcher_${download_version} \\
            $args2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher.py --version 2>&1)| sed 's/fusioncatcher.py //')
    END_VERSIONS
    """
}
