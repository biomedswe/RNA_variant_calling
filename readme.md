# RNA Variant Calling Pipeline

This project provides a **Nextflow-based pipeline** for performing variant calling on RNA-seq FASTQ files. It is designed to be modular, reproducible, and scalable for high-throughput analysis on Unix-based systems.

## 📦 Features

- STAR alignment of RNA-seq reads
- Variant calling using GATK best practices
- Support for paired-end and single-end FASTQ files
- Conda-based reproducibility with `environment.yml`

## 📁 Project Structure

Scripts/ # Nextflow scripts and pipeline logic
Data/ # Input files like sample metadata
Processed_files/ # Output directory (ignored by Git)
Logs/ # Pipeline logs and diagnostics
environment.yml # Conda environment for reproducibility


## 🧪 Requirements

- [Nextflow](https://www.nextflow.io/) >= 22.0
- [Conda](https://docs.conda.io/) (Miniconda or Anaconda)
- Access to reference genome files (e.g., STAR genome index)

## 🔁 Reproducibility

Create the Conda environment with:

```bash
conda env create -f environment.yml
conda activate variant_calling

