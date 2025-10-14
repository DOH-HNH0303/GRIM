<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-grim_logo_dark.png">
    <img alt="nf-core/grim" src="docs/images/nf-core-grim_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/grim/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/grim/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/grim/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/grim/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/grim/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/grim)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23amrlocator-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/amrlocator)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/grim** (Genomic Resistance Identification and Mapping) is a bioinformatics pipeline that maps Phoenix AMR genes onto ONT complete genomes to determine chromosomal vs plasmid location. It leverages existing Phoenix output files (GAMMA and AMRFinder results) and maps the identified AMR genes from Phoenix assemblies (derived from Illumina data) onto high-quality ONT complete genomes to provide accurate chromosomal vs plasmid classification.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Parse Phoenix AMR output files (GAMMA and AMRFinder results)
2. Extract AMR gene sequences from Phoenix assemblies (Illumina-derived)
3. Map AMR gene sequences to ONT complete genomes using BLAST
4. Classify AMR gene locations on ONT genomes (chromosomal vs plasmid)
5. Generate comprehensive AMR location mapping reports
6. Present QC and results summary ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Samplesheet Input

The pipeline supports two input formats for maximum flexibility:

#### Format 1: Individual File Paths

Specify individual paths to each Phoenix output file:

`samplesheet.csv`:

```csv
sample,gamma_ar_file,amrfinder_report,phoenix_assembly_fasta,ont_complete_genome
SAMPLE_1,/path/to/phoenix/SAMPLE_1/SAMPLE_1_ResGANNCBI_20230417_srst2.gamma,/path/to/phoenix/SAMPLE_1/SAMPLE_1_amrfinder_report.tsv,/path/to/phoenix/SAMPLE_1/SAMPLE_1_scaffolds.fasta,/path/to/ont/SAMPLE_1/SAMPLE_1_complete_genome.fasta
SAMPLE_2,/path/to/phoenix/SAMPLE_2/SAMPLE_2_ResGANNCBI_20230417_srst2.gamma,/path/to/phoenix/SAMPLE_2/SAMPLE_2_amrfinder_report.tsv,/path/to/phoenix/SAMPLE_2/SAMPLE_2_scaffolds.fasta,/path/to/ont/SAMPLE_2/SAMPLE_2_complete_genome.fasta
```

#### Format 2: Phoenix Output Directory (Recommended)

Specify the Phoenix output directory and let GRIM automatically find the required files:

`samplesheet.csv`:

```csv
sample,phoenix_outdir,ont_complete_genome
SAMPLE_1,/path/to/phoenix/output/SAMPLE_1,/path/to/ont/SAMPLE_1/SAMPLE_1_complete_genome.fasta
SAMPLE_2,/path/to/phoenix/output/SAMPLE_2,/path/to/ont/SAMPLE_2/SAMPLE_2_complete_genome.fasta
```

> [!TIP]
> Format 2 is recommended as it's simpler and automatically handles Phoenix file naming conventions. The `phoenix_outdir` should be the directory that was specified as `--outdir` when running Phoenix.

### Column Descriptions

- **`sample`**: Unique sample identifier (no spaces allowed)
- **`phoenix_outdir`**: Path to Phoenix output directory for this sample
- **`gamma_ar_file`**: Path to Phoenix GAMMA antimicrobial resistance file (`.gamma`)
- **`amrfinder_report`**: Path to Phoenix AMRFinder report (`.tsv` or `.txt`)
- **`phoenix_assembly_fasta`**: Path to Phoenix assembly file (`.fasta`, `.fa`, or `.fna`)
- **`ont_complete_genome`**: Path to ONT complete genome file (`.fasta`, `.fa`, or `.fna`, optionally gzipped)

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/grim \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/amrlocator/usage) and the [parameter documentation](https://nf-co.re/amrlocator/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/amrlocator/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/amrlocator/output).

## Credits

nf-core/grim was originally written by Seqera AI.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#grim` channel](https://nfcore.slack.com/channels/grim) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/amrlocator for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
