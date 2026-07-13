# Roche-AXELIOS/xoosnf

[![GitHub Actions CI Status](https://github.com/Roche-AXELIOS/xoosnf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/Roche-AXELIOS/xoosnf/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/Roche-AXELIOS/xoosnf/actions/workflows/linting.yml/badge.svg)](https://github.com/Roche-AXELIOS/xoosnf/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/Roche-AXELIOS/xoosnf)

## Introduction

**Roche-AXELIOS/xoosnf** is an end-to-end Nextflow pipeline for processing sequencing data produced by the Roche AXELIOS 1 sequencer.

Please see our complete [documentation](https://roche-axelios.gitbook.io/xoos/analysis-pipeline/xoosnf) for more information.

## Usage

See the [documentation](https://roche-axelios.gitbook.io/xoos/analysis-pipeline/xoosnf) for full usage instructions.

First, install all required dependencies:

- [Java >=17](https://adoptopenjdk.net/releases.html)
- [docker](https://docs.docker.com/engine/install/)

Next, download the test datasets and nextflow:

```bash
curl -sSL https://github.com/Roche-AXELIOS/xoosnf-test-datasets/releases/download/1.1.0/xoosnf-test-datasets.tar.gz | tar -xzf -

curl -sSL -o nextflow https://github.com/nextflow-io/nextflow/releases/download/v26.04.6/nextflow-26.04.6-dist
chmod +x nextflow
```

Now, you can run the pipeline from basecall FASTQs, demultiplexed FASTQs, or BAMs using:

```bash
NXF_FILE_ROOT=$(pwd)/xoosnf-test-datasets \
  ./nextflow run Roche-AXELIOS/xoosnf \
    -profile docker,test,germline_wgs_duplex \
    -params-file xoosnf-test-datasets/resources/resources.json \
    --input basecall_fastq/runsheet.csv \
    --outdir $(pwd)/output/basecall_fastq \
    --aligner bwa \
    --enable_pangenome_personalization false \
    --use_gatk_haplotype_caller true \
    --enable_variant_calling_partitioning false \
    --enable_contamination_estimation false

NXF_FILE_ROOT=$(pwd)/xoosnf-test-datasets \
  ./nextflow run Roche-AXELIOS/xoosnf \
    -profile docker,test,germline_wgs_duplex \
    -params-file xoosnf-test-datasets/resources/resources.json \
    --input bam/runsheet.csv \
    --outdir $(pwd)/output/demultiplexed \
    --aligner bwa \
    --enable_pangenome_personalization false \
    --use_gatk_haplotype_caller true \
    --enable_variant_calling_partitioning false \
    --enable_contamination_estimation false

NXF_FILE_ROOT=$(pwd)/xoosnf-test-datasets \
  ./nextflow run Roche-AXELIOS/xoosnf \
    -profile docker,test,germline_wgs_duplex \
    -params-file xoosnf-test-datasets/resources/resources.json \
    --input bam/runsheet.csv \
    --outdir output/bam \
    --aligner bwa \
    --enable_pangenome_personalization false \
    --use_gatk_haplotype_caller true \
    --enable_variant_calling_partitioning false \
    --enable_contamination_estimation false
```

## Credits

We thank the following people for their extensive assistance in the development of this pipeline:

Please see [CITATION.cff](CITATION.cff).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use Roche-AXELIOS/xoosnf for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
