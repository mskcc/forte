# mskcc/forte

[![GitHub Actions CI Status](https://github.com/mskcc/forte/actions/workflows/ci.yml/badge.svg)](https://github.com/mskcc/forte/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/mskcc/forte/actions/workflows/linting.yml/badge.svg)](https://github.com/mskcc/forte/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/mskcc/forte)

## Introduction

**mskcc/forte** is a best-practice analysis pipeline for bulk RNAseq.

- **F**unctional
- **O**bservation of
- **R**NA
- **T**ranscriptome
- **E**lements/**E**xpression

### Features

1. Read pre-processing
   1. Trimming
   2. UMI extraction and deduplication
2. Alignment
3. Transcript quantification
4. Fusion calling and annotation
5. FASTQ and BAM QC
6. Fillouts

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

```bash
nextflow run /path/to/clonedrepo/main.nf \
  --input samplesheet.csv \
  --outdir <OUTDIR> \
  --genome GRCh37 \
  -profile singularity
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](docs/usage.md).

## Credits

mskcc/forte was originally written by Anne Marie Noronha <noronhaa@mskcc.org>.

We thank the following people for their extensive assistance in the development of this pipeline:

- Allison Richards <richara4@mskcc.org>
- Alexandria Pinto <pintoa1@mskcc.org>
- Yixiao Gong <gongy@mskcc.org>

We also thank the following contributors:

- Sam Tischfield <tischfis@mskcc.org>
- Martina Bradic <bradicm@mskcc.org>
- Jun Woo <wooh@mskcc.org>
- Mark Donoghue <donoghum@mskcc.org>
- Helen Won

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use mskcc/forte for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
