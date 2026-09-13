# RNA Variant Calling Pipeline

**Created 2025-06-13 by Jonas Andersson**

This project provides a **Nextflow-based pipeline** for variant calling and allele-specific expression analysis from RNA-seq FASTQ files. It is designed to be modular, reproducible, and scalable for high-throughput analysis on Unix-based systems.

## 📦 Features

* Concatenation of FASTQ files from multiple sequencing lanes
* Quality control using FastQC and MultiQC
* Adapter and quality trimming using Trim Galore
* rRNA removal using SortMeRNA
* STAR alignment of RNA-seq reads
* RNA-seq variant calling using GATK
* WASP-based correction of allele mapping bias
* Allele-specific expression analysis using GATK ASEReadCounter
* Gene-level read quantification using featureCounts
* Conda-based reproducibility with `environment.yml`

## 📁 Project Structure

```text
Main.nf              # Main Nextflow workflow
nextflow.config      # Pipeline configuration and parameters
modules/             # Nextflow modules and helper functions
Data/                # Input files and reference data
Processed_files/     # Pipeline output (ignored by Git)
Logs/                # Pipeline logs and diagnostics
environment.yml      # Conda environment for reproducibility
```

## 🧪 Requirements

* [Nextflow](https://www.nextflow.io/) >= 22.0
* [Conda](https://docs.conda.io/) (Miniconda or Anaconda)
* [GENCODE GRCh38 release 43](https://www.gencodegenes.org/human/release_43.html)
* [SortMeRNA default rRNA reference](https://github.com/sortmerna/sortmerna/releases/download/v7.0.0/smr_v4.3_default_db.fasta.gz)
* [WASP](https://github.com/bmvdgeijn/WASP)

Clone WASP using:

```bash
git clone https://github.com/bmvdgeijn/WASP.git
```

Additional reference files required by GATK, including known variant sites, must also be downloaded and specified in the pipeline configuration.

## 🔁 Reproducibility

Clone this repository and enter the project directory:

```bash
git clone <repository-url>
cd <repository-name>
```

Create the Conda environment:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate variant_calling
```

The exact environment name is defined by the `name:` field in `environment.yml`.

## ⚙️ Configuration

Before running the pipeline, edit `nextflow.config` and specify the paths to the required input data and reference files, including:

* Input FASTQ directory
* GRCh38 reference genome
* GENCODE GTF annotation
* SortMeRNA rRNA reference
* SortMeRNA index directory
* GATK known-sites resources
* WASP installation directory
* Output directories

## ▶️ Running the Pipeline

Run the pipeline with:

```bash
nextflow run Main.nf
```

Nextflow will execute the workflow using the parameters defined in `nextflow.config`.
