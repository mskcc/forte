process PREPARE_RRNA {

    input:
    path(gtf)

    output:
    path "rna.bed", emit: rRNA_bed

    script:
    """
    ${"${gtf}".endsWith(".gz") ? "z" : ""}grep "rRNA" ${gtf} |\\
        awk '\$3 == "transcript"' | \\
        cut -f1,4,5,7,9 | \\
        perl -lane '
            /transcript_id "([^"]+)"/ or die "no transcript_id on \$.";
            print join "\t", (@F[0,1,2,3], \$1)
        ' | \\
        sort -k1V -k2n -k3n \\
        >> rna.bed

    """
}
