# FORTE

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
5. QC aggregation and reporting

## Usage

:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

Now, you can run the pipeline using:

```bash
git clone git@github.com:mskcc/forte.git
nextflow run forte/main.nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

For more details and further functionality, please refer to the [usage documentation](docs/usage.md).

## Pipeline output

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

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
<!-- If you use  mskcc/forte for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
