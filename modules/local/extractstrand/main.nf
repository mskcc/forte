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
    tuple val(meta), path("*.strandedness.tsv"), emit: strand
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/local/bin/python

    import pandas as pd

    df = pd.read_csv("${metrics}", skiprows=(lambda x: x not in [6, 7]), sep="\\t")
    df = df.drop(columns = [col for col in list(df) if col not in ['NUM_R1_TRANSCRIPT_STRAND_READS','NUM_R2_TRANSCRIPT_STRAND_READS']])

    df['input_strandedness']    = "${meta.auto_strandedness ? "auto" : meta.strandedness}"
    df['inferred_strandedness'] = df.apply(lambda row: "yes" if row['NUM_R1_TRANSCRIPT_STRAND_READS']/3 >= row['NUM_R2_TRANSCRIPT_STRAND_READS'] else "reverse" if row['NUM_R2_TRANSCRIPT_STRAND_READS']/3 >= row['NUM_R1_TRANSCRIPT_STRAND_READS'] else "no", axis=1)
    df['input_strand_correct']  = df.apply(lambda row: True if "${meta.strandedness}" == row["inferred_strandedness"] else False, axis=1)
    df.index                    = ['${meta.id}']

    desired_column_order = ['input_strandedness', 'inferred_strandedness', 'input_strand_correct','NUM_R1_TRANSCRIPT_STRAND_READS','NUM_R2_TRANSCRIPT_STRAND_READS']
    df = df[desired_column_order]

    df.to_csv("${prefix}.strandedness.tsv",sep="\\t", index=True)

    with open("versions.yml", 'w') as f:
        f.write("${task.process}:\\n")
        f.write("    pandas: 1.5.2\\n")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/python

    with open("${prefix}.strandedness.txt", 'w') as f:
        pass

    with open("versions.yml", 'w') as f:
        f.write("${task.process}:\\n")
        f.write("    pandas:1.5.2\\n")
    """
}
