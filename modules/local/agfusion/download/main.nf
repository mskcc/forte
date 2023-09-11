process AGFUSION_DOWNLOAD {
    label 'process_low'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda 'bioconda::agfusion=1.252'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/agfusion:0.0.6' :
        'cmopipeline/agfusion:0.0.6' }"

    input:
    val(ensembl_release)
    val(genome)

    output:
    path "agfusion.*.db"  , emit: agfusion_db
    path "pyensembl_cache", emit: pyensembl_cache
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def agfusion_genome = ['GRCh37','smallGRCh37','hg19'].contains(genome) ? 'hg19' :
        ['GRCh38','hg38'].contains(genome) ? 'hg38' :
        ['GRCm38','mm10'].contains(genome) ? 'mm10' : ''
    def pyensembl_species = ['GRCm38','mm10'].contains(genome) ? 'mus_musculus' : 'homo_sapiens'
    if (ensembl_release < 93) {
        """
        export PYENSEMBL_CACHE_DIR=\$PWD/pyensembl_cache

        pyensembl install --species ${pyensembl_species} --release ${ensembl_release}

        agfusion download -g ${agfusion_genome}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            agfusion: \$(agfusion -v) (fork)
        END_VERSIONS
        """
    } else {
        """
        export PYENSEMBL_CACHE_DIR=\$PWD/pyensembl_cache

        pyensembl install --species ${pyensembl_species} --release ${ensembl_release}

        curl http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/pfamA.txt.gz > pfamA.txt.gz
        gunzip pfamA.txt.gz
        agfusion build --dir . --species ${agfusion_genome} --release ${ensembl_release} --pfam pfamA.txt
        rm pfamA.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            agfusion: \$(agfusion -v) (fork)
        END_VERSIONS
        """
    }
}
