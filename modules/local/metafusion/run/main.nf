process METAFUSION_RUN {
    tag "$meta.id"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/metafusion:0.0.8' :
        'docker.io/cmopipeline/metafusion:0.0.8' }"

    input:
    tuple val(meta), path(cff)
    path genebed
    path info
    path fasta
    path blocklist
    path transcript_allowlist

    output:
    tuple val(meta), path("*final*cluster")             , emit: cluster
    tuple val(meta), path("*.exons")                    , emit: cff
    tuple val(meta), path("cis-sage.cluster")           , emit: cis
    tuple val(meta), path("problematic_chromosomes.cff"), emit: problem_chrom
    tuple val(meta), path("filters.txt")                , emit: filters
    // tuple val(meta), path("*")                       , emit: all
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ""
    def sample = "${meta.sample}"
    """
    if [ -s $cff ]; then
        export TMPDIR=\$TMPDIR
        Metafusion_forte.sh \\
            --cff $cff \\
            --outdir . \\
            --gene_bed $genebed \\
            --gene_info $info \\
            --genome_fasta $fasta \\
            --recurrent_bedpe $blocklist \\
            --clinical_genes $transcript_allowlist \\
            ${args}
    else
        echo "No fusions found by callers, returning empty files"
        touch filters.txt
        touch problematic_chromosomes.cff
        touch cis-sage.cluster
        touch empty.exons
        touch final.empty.cluster
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Metafusion: \$METAFUSION_TAG
        Metafusion_forte.sh: 0.0.2
    END_VERSIONS
    """
}
