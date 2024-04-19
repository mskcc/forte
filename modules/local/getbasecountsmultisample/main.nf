process GETBASECOUNTSMULTISAMPLE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bamtools=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/getbasecountsmultisample:1.2.4' :
        'docker://cmopipeline/getbasecountsmultisample:1.2.4' }"

    input:
    tuple val(meta), path(bam), path(bai), path(variants)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*.gbcms.maf"), optional: true, emit: maf
    tuple val(meta), path("*.gbcms.vcf"), optional: true, emit: vcf
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--omaf") ? "maf" : "vcf"
    def variant_input = variants.getExtension() == "maf" ? "--maf ${variants}" : "--vcf ${variants}"
    """
    GetBaseCountsMultiSample \\
        --thread ${task.cpus} \\
        --fasta ${fasta} \\
        ${variant_input} \\
        --bam ${prefix}:${bam} \\
        --output ${variants.getBaseName()}.gbcms.${extension} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GetBaseCountsMultiSample: \$(GetBaseCountsMultiSample -h | grep ^GetBaseCountsMultiSample | cut -f 2 -d" " )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gbcms.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GetBaseCountsMultiSample: \$(GetBaseCountsMultiSample -h | grep ^GetBaseCountsMultiSample | cut -f 2 -d" " )
    END_VERSIONS
    """
}
