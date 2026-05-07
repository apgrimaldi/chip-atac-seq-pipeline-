# Chip-atac-seq-pipeline

-----

# ATAC & ChIP-seq Analysis Pipeline

**Nextflow DSL2 pipeline for automated analysis of ChIP-seq and ATAC-seq data.**
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)


[](https://www.nextflow.io/)
[](https://www.docker.com/)

## Introduction

This pipeline is designed to process chromatin sequencing data starting from raw files (`FASTQ`) through to peak calling and annotation.

The workflow is extremely flexible: it automatically detects whether samples are **Single-End** or **Paired-End** based on the samplesheet content and adjusts MACS3 parameters accordingly.

## Usage

The pipeline can be executed directly from GitHub. Nextflow will automatically handle code download and container management.

```bash
nextflow run apgrimaldi/chip-atac-seq-pipeline- \
    -latest \
    -profile docker \
    --input samplesheet.csv \
    --protocol chip \
    --genome GRCh38 \
    --outdir "results"
```

### Main Parameters:

  * `-latest`: Forces the download of the latest code version from GitHub.
  * `-profile docker`: Runs every tool within a dedicated container (recommended).
  * `--protocol`: Defines the type of analysis (`chip` or `atac`).
  * `--genome`: Specifies the reference genome (e.g., `GRCh38`,`hg38`, ...).
  * `--input`: Path to the samplesheet CSV file.

### Genome & Annotation (Custom Genome)
 
* `--fasta_file` | Path to the reference genome FASTA file (e.g., `.fna`, `.fasta`). 
* `--gtf_file` | Path to the gene annotation file (GTF format) for QC and peak annotation. 
* `--macs_gsize` | Effective genome size for MACS3 (e.g., `2.7e9` for human, `hs`, or `mm`). 
* `--blacklist` | Path to the BED file containing blacklisted regions to be excluded. 

## Usage Example for Custom Genome

To run the pipeline with a custom genome and all required annotations, use the following command:

```bash
nextflow run main.nf \
  -profile docker \
  --protocol <chip/atac> \
  --input path/to/samplesheet.csv \
  --fasta_file "/path/to/reference_genome.fna" \
  --gtf_file "/path/to/annotation.gtf" \
  --macs_gsize <genome_size> \
  --blacklist "/path/to/blacklist.bed" \
  -resume


## Pipeline Summary

The workflow performs the following steps:

1.  **Quality Control**: Raw read quality control using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2.  **Trimming**: Removal of adapters and low-quality bases with [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
3.  **Alignment**: Read mapping to the reference genome via [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
4.  **Alignment Management**: Processing, sorting, and indexing BAM files using [SAMtools](http://www.htslib.org/).
5.  **Duplicates Management**: Identification and removal of PCR duplicates with [Picard MarkDuplicates](https://broadinstitute.github.io/picard/).
6.  **Blacklist Removal**: Filtering out reads overlapping problematic genomic regions using [BEDTools](https://bedtools.readthedocs.io/).
7.  **Quality Metrics & Statistics**: Generating comprehensive alignment statistics with [SAMtools stats](http://www.htslib.org/doc/samtools-stats.html).
8.  **Enrichment Analysis**: Generating Fingerprint to assess IP strength using [deepTools](https://deeptools.readthedocs.io/).
9.  **Peak Calling**: Identification of enriched genomic regions (Narrow/Broad) with [MACS3](https://github.com/macs3-project/MACS).
10. **Annotation**: Functional annotation of peaks relative to gene features using [HOMER](http://homer.ucsd.edu/homer/).
11. **QC Metrics**: Calculation of the Fraction of Reads in Peaks (FRiP score).
12. **MultiQC**: Compilation of an interactive report aggregating stats from all steps using [MultiQC](https://multiqc.info/).

-----

## Input (Samplesheet)

The `samplesheet.csv` file must be formatted as follows:

```csv
sample,fastq_1,fastq_2,antibody,control
IP_gH2AX_DOXO_1_S19_R1_001,data/IP_gH2AX_DOXO_1_S19_R1_001.fastq.gz,,IgG,IP_IgG_DOXO_1_S22_R1_001
IP_IgG_DOXO_1_S22_R1_001,data/IP_IgG_DOXO_1_S22_R1_001.fastq.gz,,,
```

The columns must be structured as follows:

  * **sample**: Unique name for the sample.
  * **fastq\_1**: Full path to FastQ file 1 (Read 1). Must end in `.fastq.gz` or `.fq.gz`.
  * **fastq\_2**: Full path to FastQ file 2 (Read 2). Leave **empty** for Single-End samples.
  * **antibody**: Name of the antibody used (e.g., `H3K27me3`).
  * **control**: Name of the sample to be used as control (must match a value in the `sample` column).

## Output

Results are organized in the `results/` folder:

  * **`00_MultiQC/`**: Interactive HTML report (including genome info and tool versions).
  * **`01_fastqc/`**: Initial read quality.
  * **`04_alignment/`**: Filtered and indexed BAM files.
  * **`05_peaks/`**: `.narrowPeak` / `.broadPeak` files and MACS3 logs.
  * **`06_bigwig/`**: `.bw` files for track visualization.
  * **`07_homer_annotation/`**: Gene annotation tables.

https://github.com/apgrimaldi/chip-atac-seq-pipeline-s

## Credits

Developed with passion by **Annapaola** (@apgrimaldi).

> *Note: The MultiQC report is configured to display configuration information (Genome/Protocol) at the top and the software versions summary at the bottom, following nf-core standards.*

-----
