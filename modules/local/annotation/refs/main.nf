process REFS_ANNOTATIONS {
	label "process_single"

	input:
	val(genome)
	val(ensembl_release)
	path(target_genes)
	path(gtf)
	path(ucsc_pfam)

	output:
	path "reference_transcripts.txt"			,	emit: reference_transcripts
	path "svtable.txt"							, 	emit: sv_table
	path "refFlat.txt"							,	emit: refFlat
	path "refFlat_summary.txt"					,	emit: refFlat_summary
	path "Pfam_kinase_domains.txt"				,	emit: kinase_domains

	when:
	task.ext.when == null || task.ext.when

	script:
	"""
	prep_annot_refs.sh \\
		-t $target_genes \\
		-g $gtf \\
		-p $ucsc_pfam \\
		-r reference_transcripts.txt \\
		-s svtable.txt \\
		-f refFlat.txt \\
		-u refFlat_summary.txt \\
		-k Pfam_kinase_domains.txt

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
		python: 3.12
	END_VERSIONS
	"""
}