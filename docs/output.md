# mskcc/forte: Output

## Introduction

This document describes the output produced by the FORTE pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Read Preprocessing](#read-preprocessing)
- [Alignment](#alignment)
- [Quantification](#quantification)
- [Fusion Calling](#fusion-calling)
- [Fusion Merging and Annotation](#fusion-merging-and-annotation)
- [QC](#qc)
- [Fillouts](#fillouts)
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Read Preprocessing

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/fastp/`
  - `*.fastp.html`
  - `*.fastp.json`
  - `*.fastp.log`
  - `*.fastp.fastq.gz`
- `analysis/<sample>/umitools/extract/`
  - `logs/.umi_extract.log`

</details>

[FastP](https://github.com/OpenGene/fastp) gives general quality metrics about your sequenced reads and also trims the reads according to base quality and presence of adapter sequences.

[UMI-tools extract](https://umi-tools.readthedocs.io/en/latest/reference/extract.html) removes UMI sequences from reads and adds it to the read header. As a result, aligners do not attempt to align the UMI sequence and the aligned reads will be ready for deduplication.

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/STAR/`
  - `*.Aligned.sortedByCoord.out.bam`
  - `*.Aligned.sortedByCoord.out.bam.bai`
  - `log/`
    - `*.Log.out`
    - `*.Log.final.out`
    - `*.Log.progress.out`
    - `*.ReadsPerGene.out.tab`
    - `*.SJ.out.tab`
- `analysis/<sample>/umitools/dedup/`
  - `*.dedup.bam`
  - `*.dedup.bam.bai`
  - `logs/`
    - `*.dedup_edit_distance.tsv`
    - `*.dedup_per_umi_per_position.tsv`
    - `*.dedup_per_umi.tsv`

</details>

[STAR](https://github.com/alexdobin/STAR) is an ultrafast universal RNA-seq aligner.

[UMI-tools dedup](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) deduplicates reads based on the mapping co-ordinate and the UMI attached to the read.

### Quantification

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/featurecounts/`
  - `*.gene.featureCounts.txt`
- `analysis/<sample>/kallisto/`
  - `abundance.h5`
  - `abundance.tsv`
  - `run_info.json`
  - `*.log.txt`

</details>

[featureCounts](https://subread.sourceforge.net/featureCounts.html) takes a file with aligned sequencing reads, plus a list of genomic features and counts how many reads map to each feature.

[Kallisto](http://pachterlab.github.io/kallisto/) quantifies abundances of transcripts from RNA-Seq data using high-throughput sequencing reads.

### Fusion Calling

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/arriba/`
  - `*.fusions.discarded.tsv`
  - `*.fusions.tsv`
  - `*_arriba.cff`
- `analysis/<sample>/fusioncatcher/`
  - `*.fusioncatcher.fusion-genes.hg19.txt`
  - `*.fusioncatcher.fusion-genes.txt`
  - `*.fusioncatcher.log`
  - `*.fusioncatcher.summary.txt`
  - `*_fusioncatcher.cff`
- `analysis/<sample>/starfusion/`
  - `*.starfusion.abridged.coding_effect.tsv`
  - `*.starfusion.abridged.tsv`
  - `*.starfusion.fusion_predictions.tsv`
  - `*_starfusion.cff`
  - `STAR/`
    - `*.Chimeric.out.junction`
    - `log/`
      - `*.Log.final.out`
      - `*.Log.out`
      - `*.Log.progress.out`
      - `*.SJ.out.tab`

</details>

[Arriba](https://arriba.readthedocs.io/en/latest/) uses the STAR aligner to detect of gene fusions from RNA-Seq data.

[FusionCatcher](https://github.com/ndaniel/fusioncatcher) searches for novel/known somatic fusion genes, translocations, and chimeras in RNA-seq data.

[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) uses the STAR aligner to identify candidate fusion transcripts supported by Illumina reads.

[CommonFusionFormat (CFF)](https://github.com/ccmbioinfo/MetaFusion/wiki/metafusion-file-formats) file is originally created by the developers of the [MetaFusion](https://github.com/mskcc/MetaFusion) tool. Forte generates `*.cff` files that can be used with MetaFusion.

### Fusion Merging and Annotation

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/metafusion`
  - `*.final.cff`
  - `*.unfiltered.cff`
  - `intermediates/`
    - `cis-sage.cluster`
    - `*.cff.cleaned_chr.renamed.reann.WITH_SEQ.exons`
    - `*_metafusion_cluster.unfiltered.cff`
    - `final.n1.cluster`
    - `problematic_chromosomes.cff`

</details>

FORTE uses a custom fork of [Metafusion](https://github.com/mskcc/MetaFusion) to filter, cluster and annotate the fusion calls. Several `intermediate` files are included in the output.

`Fusion_effect` information is added using a custom fork of [AGFusion](https://github.com/anoronh4/AGFusion).

[FusionAnnotator.py from the oncokb-annotator](https://github.com/oncokb/oncokb-annotator/blob/master/FusionAnnotator.py) is also run and added to the final cff file.

### QC

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/picard/`
  - `*.rna_metrics`
  - `*.CollectHsMetrics.coverage_metrics`
- `analysis/<sample>/rseqc/`
  - `*.bam_stat.txt`
  - `*.DupRate_plot.pdf`
  - `*.DupRate_plot.r`
  - `*.infer_experiment.txt`
  - `*.inner_distance_freq.txt`
  - `*.inner_distance_mean.txt`
  - `*.inner_distance_plot.pdf`
  - `*.inner_distance_plot.r`
  - `*.inner_distance.txt`
  - `*.junction_annotation.log`
  - `*.junction.bed`
  - `*.junction.Interact.bed`
  - `*.junction_plot.r`
  - `*.junctionSaturation_plot.pdf`
  - `*.junctionSaturation_plot.r`
  - `*.junction.xls`
  - `*.pos.DupRate.xls`
  - `*.read_distribution.txt`
  - `*.seq.DupRate.xls`
  - `*.splice_events.pdf`
  - `*.splice_junction.pdf`
- `analysis/<sample>/multiqc/`
  - `dedupbam_multiqc_report_data/`
    - `*.json`
    - `*.log`
    - `*.txt`
  - `dedupbam_multiqc_report.html`
  - `dedupbam_multiqc_report_plots/`
    - `pdf/*.pdf`
    - `png/*.png`
    - `svg/*.svg`
  - `dupbam_multiqc_report_data/`
    - `*.json`
    - `*.log`
    - `*.txt`
  - `dupbam_multiqc_report.html`
  - `dupbam_multiqc_report_plots/`
    - `pdf/*.pdf`
    - `png/*.png`
    - `svg/*.svg`

</details>

[Picard's CollectHsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard-) collects hybrid-selection (HS) metrics for a SAM or BAM file. This is only produced if baitset is indicated in the samplesheet.

[Picard's CollectRnaSeqMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard-) produces RNA alignment metrics for a SAM or BAM file.

[RSeQC](https://rseqc.sourceforge.net/) provides a number of useful modules that can comprehensively evaluate high throughput RNAseq data.

[MultiQC](https://multiqc.info/) is a visualization tool that searches a given directory for analysis/qc logs and compiles a HTML report. Most of the pipeline QC results are visualized in the report and further statistics are available in the report data directory. FORTE produces a second MultiQC report for each sample that has UMI. FORTE also produces 1-2 reports under the `multiqc/` folder where all samples are aggregated together, one for non-deduplicated results and the other for deduplicated results.

### Fillouts

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/fillouts`
  - `*.fillout.maf`

</details>

[GetBaseCountsMultiSample (GBCMS)](https://github.com/zengzheng123/GetBaseCountsMultiSample) calculates the base counts in a given BAM file for all the sites in a given MAF file

FORTE uses a custom script to output a MAF file that combines all original columns and new columns from GBCMS.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
