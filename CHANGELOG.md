# mskcc/forte: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## dev

### `Added`

- [#117](https://github.com/mskcc/forte/pull/117) - add supporting-reads_gene-fusions\*.zip files to fusioncatcher outputs

- [#118](https://github.com/mskcc/forte/pull/118) - change the way the plug-n-play starfusion reference is downloaded.

- [#126](https://github.com/mskcc/forte/pull/126) - enable transcript prioritization in Metafusion

- [#128](https://github.com/mskcc/forte/pull/128) - full support for GRCh38 added

- [#141](https://github.com/mskcc/forte/pull/141) - Add portcullis

- [#138](https://github.com/mskcc/forte/pull/128) - enable clinical gene expansion in agfusion

### `Fixed`

- [#119](https://github.com/mskcc/forte/pull/119) - change script error behavior in METAFUSION_RUN process

- [#125](https://github.com/mskcc/forte/pull/125) - update upload-artifact version because the version previously in use (v2) is deprecated.

- [#124](https://github.com/mskcc/forte/pull/124) - ensure genebed file as 0based start site

- [#127](https://github.com/mskcc/forte/pull/127) - allow dynamic increase of memory for process_single label

- [#132](https://github.com/mskcc/forte/pull/132) - fix generate cff split/span logic for fusioncatcher and arriba

- [#133](https://github.com/mskcc/forte/pull/133) - Template update from nf-core/tools v3.1.1, including addition of pipeline initialization and completion workflows, and input schemas. Also changing base branch to dev for more harmony with the template.

- [#135](https://github.com/mskcc/forte/pull/135) - sort read_group and fastq_pair_id values before concatenation

- [#139](https://github.com/mskcc/forte/pull/139) - Add nftests and CI workflow

- [#142](https://github.com/mskcc/forte/pull/142) - Template update from nf-core/tools v3.2.0

- [#148](https://github.com/mskcc/forte/pull/148) - Consolidate arriba bam alignment and primary bam alignment into a single process

- [#149](https://github.com/mskcc/forte/pull/149) - Temporary downgrade of nf-schema to 2.2.0 to allow tests to complete

- [#151](https://github.com/mskcc/forte/pull/151) - Replace samtools bam2fq with gatk4 samtofastq for better handling of paired reads in a position sorted bam

- [#159](https://github.com/mskcc/forte/pull/159) - Exclude any fusion with NA as one of the genes from agfusion clinical run

- [#154](https://github.com/mskcc/forte/pull/151) - Fix publishing of sample-level and batch-level multiqc reports and add parameter to turn off plot export in multiqc

- [#160](https://github.com/mskcc/forte/pull/160) - Fix grouping of multiple pairs of fastqs in a single sample

### `Dependencies`

### `Deprecated`

## v1.0.0 - 2023-09-23

Initial release of mskcc/forte, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
