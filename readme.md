# RNA Variant Calling Pipeline
# Created 2025-06-13 by Jonas Andersson

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

✂️ Trim Galore Summary

Trimming mode: paired-end
Trim Galore version: 0.6.10
Cutadapt version: 2.6
Python version: 3.7.12          
Number of cores used for trimming: 8
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp

Note: Trim Galore was executed outside the main Conda environment due to dependency conflicts.
