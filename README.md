# ![nf-core/rrms](docs/images/nf-core-rrms_logo_light.png#gh-light-mode-only) ![nf-core/rrms](docs/images/nf-core-rrms_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rrms/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


## Introduction

SkewX is a nextflow pipeline to measure skewed X inactivation from long-read sequencing of native DNA, either from Pacbio or Nanopore.
It starts from bam files that include modified basecalls for 5mCG. It first calls heterozygous variants with DeepVariant and phases them into haplotypes with WhatsHap. Then it also clusters reads based on their methylation profile over CpG islands.

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-core/rrms** is a bioinformatics best-practice analysis pipeline for qg.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/rrms/results).

## Pipeline summary

<!-- Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Basecall raw nanopore data (.pod5) with (['Dorado'](https://github.com/nanoporetech/dorado))
2. Call variants with (['DeepVariant'](https://github.com/google/deepvariant))
3. Call structural variants with (['Sniffles2'](https://github.com/fritzsedlazeck/Sniffles))
3. Phase SNPs with (['WhatsHap'](https://whatshap.readthedocs.io/en/latest/index.html))
4. Haplotype and tag reads with (['WhatsHap'](https://whatshap.readthedocs.io/en/latest/index.html))
5. Create bigwigs for easy methylation visualisation in IGV with (['modbam2bed'](https://github.com/epi2me-labs/modbam2bed));(['bedGraphToBigWig'](https://genome.ucsc.edu/goldenPath/help/bigWig.html))
6. Create ROI methylation plots with (['NanoMethViz'](https://www.bioconductor.org/packages/release/bioc/html/NanoMethViz.html))

## Quick Start

1. Install or module load [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. IMPORTANT - ensure you mount singularity to your home directory (include "export NXF_SINGULARITY_HOME_MOUNT=true" in your .bashrc or to your session environment before launching pipeline - by default Singularity will not be able to find your home)

4. Ensure required files (.bed files, .fa reference) are properly specified as parameters in the config (nextflow.config)

5. For now, make sure modbam2bed and bedGraphToBigWig binaries are cloned and made in projectDir (TODO: create process for their installation)

6. Download the pipeline and test it on a minimal dataset (3-4 .pod5s) with a single command:

   ```bash
   nextflow main.nf --input sample_sheet.csv --outdir <OUTDIR> --fasta <REFERENCE.fa> -profile singularity
   ```
7. Add any additional regions of interest for NanoMethViz to plot (in the NanoMethViz module)

8. Start running your own analysis!

   ```bash
   nextflow main.nf --input samplesheet.csv --outdir rrms_feb_2024 --fasta /home/reference/chm13v2.0.fa -profile singularity
   ```

## Documentation

The nf-core/rrms pipeline comes with documentation about the pipeline [usage](https://nf-co.re/rrms/usage), [parameters](https://nf-co.re/rrms/parameters) and [output](https://nf-co.re/rrms/output).
(links broken)

## Credits

nf-core/rrms was originally written by Quentin Gouil.

We thank the following people for their extensive assistance in the development of this pipeline:

- James Lancaster (updated pipeline tools and functionality)

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rrms` channel](https://nfcore.slack.com/channels/rrms) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/rrms for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
