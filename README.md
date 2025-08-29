# SkewX


  [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.3-brightgreen.svg)](https://www.nextflow.io/)
  [![DOI]()](https://doi.org/10.1101/2024.03.20.585856)

## Introduction

**SkewX** is a nextflow pipeline to measure skewed X inactivation from long-read sequencing of native DNA, either with Pacbio or Nanopore or technologies.
It starts from bam files that include modified basecalls for 5mCG. It first calls heterozygous variants with DeepVariant and phases them into haplotypes with WhatsHap. Then it also clusters reads based on their methylation profile over CpG islands, and pools this haplotype and epiallele information to measure the skew in X inactivation.


The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

The required input is modbam files with 5mCG information. Then:

1. If the reads are not already aligned, align to the reference genome with ['Minimap2'](https://github.com/lh3/minimap2)
2. If multiple samples per individual are present, for instance multiple tissues, merge them into a single bam file
3. Call variants with ['DeepVariant'](https://github.com/google/deepvariant)
4. Phase SNPs with ['WhatsHap'](https://whatshap.readthedocs.io/en/latest/index.html)
5. Haplotype and tag reads with ['WhatsHap'](https://whatshap.readthedocs.io/en/latest/index.html)
6. Cluster reads based on methylation profile with ['NanoMethViz'](https://www.bioconductor.org/packages/release/bioc/html/NanoMethViz.html)
7. Measure skew in X inactivation and generate a report for each individual.

## Quick Start

1. Install or module load [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. IMPORTANT - ensure you mount singularity to your home directory (include "export NXF_SINGULARITY_HOME_MOUNT=true" in your .bashrc or to your session environment before launching pipeline - by default Singularity will not be able to find your home)

4. Ensure required files (.bed files, .fa reference) are properly specified as parameters in the config (nextflow.config)

5. Start running your own analysis!

   ```bash
   nextflow main.nf --input samplesheet.csv --outdir skew_results --fasta chm13v2.0.fa --cgi CGIs_CHM13v2_chrX.bed -profile singularity
   ```

<!-- TO DO: accept phased vcf to skip deepvariant, or accept haplotyped bam to skip deepvariant+whatshap? How to handle that without needing another samplesheet? I could provide haplotyped bams-->


## Documentation

## Example data

An example dataset is available in the `test_data` directory of this repository. The dataset contains a small region of the mouse X chromosome, with a BAM file with methylation information. The pipeline can be run on this dataset with the following command:

```bash
nextflow main.nf --input test_data/test_data_samplesheet.csv --outdir skew_test_results --reference test_data/mm10_chrX.fa --cgi test_data/mm10_chrX_CGI.bed -profile test
```

## Credits

SkewX was originally written by Quentin Gouil, James Lancaster and Ed Yang.

We thank the following people for their extensive assistance in the development of this pipeline:

- Kathleen Zeglinski for her superior nextflow expertise
- Shian Su for implementing new features in NanoMethViz



## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
If you use  **SkewX** for your analysis, please cite it using the following doi: [10.1101/2024.03.20.585856](https://doi.org/10.1101/2024.03.20.585856)

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
