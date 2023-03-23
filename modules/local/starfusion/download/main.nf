process STARFUSION_DOWNLOAD {
    tag 'star-fusion'

    conda "bioconda::dfam=3.3 bioconda::hmmer=3.3.2 bioconda::star-fusion=1.10.0 bioconda::trinity=date.2011_11_2 bioconda::samtools=1.9 bioconda::star=2.7.8a"
    container "docker.io/trinityctat/starfusion:1.10.1"

    input:
    val starfusion_plug_n_play_address

    output:
    path "ctat_genome_lib_build_dir/"             , emit: reference
    path "ctat_genome_lib_build_dir/ref_annot.gtf", emit: chrgtf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget ${starfusion_plug_n_play_address} --no-check-certificate

    tar xvf *.plug-n-play.tar.gz

    rm *.plug-n-play.tar.gz

    mv */ctat_genome_lib_build_dir .
    """
}
