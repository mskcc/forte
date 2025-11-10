process ANNOTATION_RUN {
	tag "$meta.id"
	label "process_single"

	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
		'docker://100813/target-fusion-annotator:1.0.1':
		'docker.io/100813/target-fusion-annotator:1.0.1' }"

	input:
	tuple val(meta), val(filtered_fusion) 
	path(ref_txs)  
	path(sv_table) 
	path(refFlat)
	path(refFlat_summary)
	path(kinase_domains)
	path(tumor_suppressor_genes)
	path(oncokb_known_fusions)
	path(keygenes)

	output:
	tuple val(meta), path("${meta.sample}_AllAnnotatedSVs.txt")	, emit: final_annotated
	tuple val(meta), path("${meta.sample}_AllAnnotatedSVs.novel.txt")	, emit: novel
	tuple val(meta), path("${meta.sample}_AllAnnotatedSVs.dropped.txt")	, emit: dropped_annotated
	path "versions.yml"

	when:
	task.ext.when == null || task.ext.when

	script:
	def args   = task.ext.args ?: ''
	def sample = "${meta.sample}"
	"""
	/app/run_fusion_annotator.sh\\
		${filtered_fusion} \\
		${sv_table} \\
		${ref_txs} \\
		${refFlat_summary} \\
		${tumor_suppressor_genes} \\
		${oncokb_known_fusions} \\
		${keygenes} \\
		${kinase_domains}

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
		python: 3.12
	END_VERSIONS
	"""
}