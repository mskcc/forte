## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**anoronh4/forte** is a best-practice analysis pipeline for bulk RNAseq.

- **F**unctional
- **O**bservation of
- **R**NA
- **T**ranscriptome
- **E**lements/**E**xpression

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Read pre-processing
   1. Trimming
   2. UMI extraction and deduplication
2. Alignment
3. Transcript quantification
4. Fusion calling and annotation
5. FASTQ and BAM QC

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run anoronh4/forte -profile test,singularity --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run anoronh4/forte --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile singularity
   ```

   The input file should contain IDs, paths and meta-data pertaining to each sample. The following is a description of each field that can be used. Fields that do not have a default value are required; those that do are not.

   | Header  | Type        | Values               | Defaults |
   | :------ | :---------- | :------------------- | :------- |
   | sample  | `str`       |                      | (none)   |
   | umi     | `str`/`int` | `NNNXX`/`3`          | `''`     |
   | umi2    | `str`/`int` | `NNNXX`/`3`          | `''`     |
   | strand  | `str`       | `yes`/`no`/`reverse` | `no`     |
   | fastq_1 | `str`       | `/path/to/*fastq.gz` | (none)   |
   | fastq_1 | `str`       | `/path/to/*fastq.gz` | (none)   |

   If you are running on juno, chain the `juno` profile (i.e. `-profile singularity,juno`) to take advantage of local resources on juno.

## Credits

anoronh4/forte was originally written by Anne Marie Noronha <noronhaa@mskcc.org>.

We thank the following people for their extensive assistance in the development of this pipeline:

- Allison Richards <richara4@mskcc.org>
- Alexandra Pinto <pintoa1@mskcc.org>
- Yixiao Gong <gongy@mskcc.org>

We also thank the following contributors:

- Sam Tischfield <tischfis@mskcc.org>
- Martina Bradic <bradicm@mskcc.org>
- Jun Woo <wooh@mskcc.org>
- Helen Won
<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
