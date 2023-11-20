process EXTRACTSTRAND {
    tag "$meta.id"
    label 'process_single'

    conda "biocontainers::pandas:1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
            'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(metrics)

    output:
    tuple val(meta), path("*.determination.txt"), emit: strand
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/local/bin/python

    import pandas as pd

    df = pd.read_csv("${metrics}", skiprows=(lambda x: x not in [6, 7]), sep="\\t")
    r1_transcript_strand_reads = int(df.iloc[0]['NUM_R1_TRANSCRIPT_STRAND_READS'])
    r2_transcript_strand_reads = int(df.iloc[0]['NUM_R2_TRANSCRIPT_STRAND_READS'])

    if r1_transcript_strand_reads/3 > r2_transcript_strand_reads:
        determination = "yes"
    elif r2_transcript_strand_reads/3 > r1_transcript_strand_reads:
        determination = "reverse"
    else:
        determination = "no"

    strandedness_correct = True
    if "${meta.strandedness}" != "auto":
        if determination == "${meta.strandedness}":
            strandedness_correct = False

    with open("${prefix}.determination.txt",'w') as f:
        f.write("${meta.id}\\t${meta.strandedness}\\t" + determination + "\\t" + str(strandedness_correct) + "\\n")

    with open("versions.yml", 'w') as f:
        f.write("${task.process}:")
        f.write("    pandas:1.5.2")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/python

    with open("${prefix}.determination.txt", 'w') as f:
        pass

    with open("versions.yml", 'w') as f:
        f.write("${task.process}:")
        f.write("    pandas:1.5.2")
    """
}
