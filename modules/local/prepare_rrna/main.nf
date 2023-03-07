process PREPARE_RRNA {

    input:
    path(gtf)
    path(refflat)

    output:
    path "rna.bed", emit: rRNA_bed

    script:
    if (gtf) {
        """
        (${"${gtf}".endsWith(".gz") ? "z" : ""}grep "rRNA" ${gtf} || true) | \\
            awk '\$3 == "transcript"' | \\
            cut -f1,4,5,7,9 | \\
            perl -lane '
                /transcript_id "([^"]+)"/ or die "no transcript_id on \$.";
                print join "\t", (@F[0,1,2,3], \$1)
            ' | \\
            (grep -vP "^HG|^HSCHR" || true) | \\
            sort -k1V -k2n -k3n \\
            > rna.bed

        """
    } else {
        """
        (${"${refflat}".endsWith(".gz") ? "z" : ""}grep -P "^RNA5|^RNA1|^RNA2" ${refflat} || true) | \\
            awk -F"\\t" -v OFS="\\t" '{ print \$3,\$5,\$6,\$4,\$2 }' | \\
            (grep -vP "^HG|^HSCHR" || true) | \\
            sort -k1V -k2n -k3n \\
            > rna.bed
        """
    }
}
