process AGFUSION_DOWNLOAD {
    label 'process_low'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda 'bioconda::agfusion=1.252'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/agfusion:0.0.7' :
        'docker.io/cmopipeline/agfusion:0.0.7' }"

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
    def pyensembl_genome  = agfusion_genome == "hg19" ? "GRCh37" :
        agfusion_genome == "hg38" ? 'GRCh38' :
        agfusion_genome == "mm10" ? 'GRCm38' :
        ''
    def pyensembl_species = ['GRCm38','mm10'].contains(genome) ? 'mus_musculus' : 'homo_sapiens'
    def pyensembl_species_capitalized = pyensembl_species[0].toUpperCase() + pyensembl_species.substring(1)
    if (ensembl_release < 112) {
        """
        export PYENSEMBL_CACHE_DIR=\$PWD/pyensembl_cache

        pyensembl install --species ${pyensembl_species} --release ${ensembl_release}

        if [ ! -f \$PYENSEMBL_CACHE_DIR/pyensembl/${pyensembl_genome}/ensembl${ensembl_release}/${pyensembl_species_capitalized}.${pyensembl_genome}.${ensembl_release}.gtf.gz ] ; then
            wget https://ftp.ensembl.org/pub/release-${ensembl_release}/gtf/${pyensembl_species}/${pyensembl_species_capitalized}.${pyensembl_genome}.${ensembl_release}.gtf.gz -O \$PYENSEMBL_CACHE_DIR/pyensembl/${pyensembl_genome}/ensembl${ensembl_release}/${pyensembl_species_capitalized}.${pyensembl_genome}.${ensembl_release}.gtf.gz
        fi

        agfusion download -s ${pyensembl_species} -r ${ensembl_release}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            agfusion: \$(agfusion -v) (fork)
        END_VERSIONS
        """
    } else {
        """
        export PYENSEMBL_CACHE_DIR=\$PWD/pyensembl_cache

        pyensembl install --species ${pyensembl_species} --release ${ensembl_release}

        if [ ! -f \$PYENSEMBL_CACHE_DIR/pyensembl/${pyensembl_genome}/ensembl${ensembl_release}/${pyensembl_species_capitalized}.${pyensembl_genome}.${ensembl_release}.gtf.gz ] ; then
            wget https://ftp.ensembl.org/pub/release-${ensembl_release}/gtf/${pyensembl_species}/${pyensembl_species_capitalized}.${pyensembl_genome}.${ensembl_release}.gtf.gz -O \$PYENSEMBL_CACHE_DIR/pyensembl/${pyensembl_genome}/ensembl${ensembl_release}/${pyensembl_species_capitalized}.${pyensembl_genome}.${ensembl_release}.gtf.gz
        fi

        curl http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/database_files/pfamA.txt.gz > pfamA.txt.gz
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
