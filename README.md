Written with the help of github-copilot (GPT-4o), this README provides detailed instructions for setting up and running the RNA-Seq processing pipeline. 

# RNA-Seq Processing Pipeline

This pipeline processes RNA-Seq data, including trimming, quality control, read mapping, and quantification. It uses `fastp` for trimming and quality control, `STAR` for read mapping, and `HTSeq` for quantification. I make no claims that this is the best piepline, but it works for my purposes and should be adaptable to your needs with minimal modifications.

---

## Table of Contents
1. [Requirements](#requirements)
2. [Setup](#setup)
3. [Pipeline Overview](#pipeline-overview)
4. [Input Files](#input-files)
5. [Output Files](#output-files)
6. [Running the Pipeline](#running-the-pipeline)
7. [Directory Structure](#directory-structure)
8. [Troubleshooting](#troubleshooting)

---

## Requirements

### Software
- **Conda/Mamba**: Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), [Anaconda](https://www.anaconda.com/products/distribution) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). Nowdays, I prefer mamba because its fast, but a modern version of conda also works. If you use mamba, just replace `conda` with `mamba` in the commands below.
- **Snakemake**: Install via Conda:
  ```bash
  conda create -n snakemake_env -c bioconda snakemake
  conda activate snakemake_env
  conda install -c bioconda snakemake
  ```

### Conda Environments
The pipeline uses the following Conda environments:
- `fastp`: For trimming and quality control.
- `qc`: For additional quality control tools (`fastqc`, `multiqc`).
- `rna_seq_processing`: For read mapping (`STAR`), BAM file processing (`samtools`), and quantification (`HTSeq`).

Environment YAML files are provided in `workflow/envs/`. See the [Setup](#setup) section for instructions on creating these environments.

---

## Setup

1. **Clone the Repository**:
   Copy the pipeline files to your working directory.

2. **Install Snakemake**:
   ```bash
   conda install -c bioconda snakemake
   ```

3. **Prepare Input Files**:
   - Place your raw FASTQ files in the directory specified in `samples.csv`.
   - Ensure the `samples.csv` file contains the correct metadata for your samples.

---

## Pipeline Overview

### Steps
1. **Trimming and Quality Control**:
   - `fastp` trims raw reads and generates quality control reports.
2. **Read Mapping**:
   - `STAR` maps trimmed reads to the reference genome.
3. **BAM File Processing**:
   - `samtools` sorts and indexes BAM files.
4. **Quantification**:
   - `HTSeq` generates count tables for downstream analysis.

### Tools
- `fastp`: Trimming and QC.
- `fastqc`/`multiqc`: Additional QC visualization.
- `STAR`: Read alignment.
- `samtools`: BAM file processing.
- `HTSeq`: Quantification.

---

## Input Files

### 1. `samples.csv`
A CSV file containing metadata for your samples. Example:

```csv
strain,biorep,lane,fq1,fq2
strain1,1,L001,data/strain1_biorep1_L001_R1.fastq.gz,data/strain1_biorep1_L001_R2.fastq.gz
strain1,1,L002,data/strain1_biorep1_L002_R1.fastq.gz,data/strain1_biorep1_L002_R2.fastq.gz
strain2,1,L001,data/strain2_biorep1_L001_R1.fastq.gz,data/strain2_biorep1_L001_R2.fastq.gz
strain2,2,L001,data/strain2_biorep2_L001_R1.fastq.gz,data/strain2_biorep2_L001_R2.fastq.gz
```
Naming convention:
- `strain`: Sample strain identifier. This will define the subfolders that results will be split into.
- `biorep`: Biological replicate number. This will be tagged onto the output file but does not affect analysis. You can also include information like time here.
- `lane`: Sequencing lane identifier. All files with matching `strain` and `biorep` will be merged together, so this will not show up in the final count files.

### 2. Reference Genome
- **Genome Index**: Directory containing the STAR genome index.
- **Annotation File**: GFF file for HTSeq.

Update the paths to these files in `config/config.yaml`.

---

## Output Files

### Key Output Directories
- **Trimmed FASTQ Files**: `results/trimmed_fq/`
- **Quality Control Reports**: `results/qc/`
- **Mapped Reads (SAM/BAM)**: `results/mapped_reads/` and `results/mapped_sorted_reads/`
- **HTSeq Counts**: `results/htseq_counts/`

---

## Running the Pipeline

1. **Edit Configuration**:
   Update `config/config.yaml` to specify the paths to your `samples.csv`, reference genome, and output directories.

2. **Edit Cluster Configuration**:
   Double check over `profile.config.yaml` to ensure setup is correct for DCC slurm environment. You can also adjust
   default resource and partition settings here.

3. **Edit slurm submission script**:
   To run with slurm, edit the provided `sample_slurm.sh` script. Note that resource allocations in this script are only for the control script managing the jobs and NOT for the jobs themselves (those are edited either in the job definitions in the Snakefile or golbally in `profile.config.yaml`). Remember to adjust the conda/mamba profile and activation lines to match your environement setup. Your first run of this script will be a dry run. Check the output file to make sure all of your jobs are being correctly identfied. Once that looks good, remove the -n flag to submit the full set of jobs to the cluster.

4. **Run Snakemake**:
   Execute the pipeline with:
   ```bash
   sbatch sample_slurm.sh
   ```
5. **Check Outputs**:
   - Trimmed FASTQ files will be in `results/trimmed_fq/`.
   - Quality control reports will be in `results/qc/`.
   - HTSeq count files will be in `results/htseq_counts/`.

---

## Directory Structure

The pipeline assumes the following directory structure:

```
project/
├── Snakefile
├── config/
│   ├── config.yaml
│   └── samples.csv
├── data/
│   ├── strain1_biorep1_L001_R1.fastq.gz
│   ├── strain1_biorep1_L001_R2.fastq.gz
│   └── ...
├── resources/
│   └── genomes/
│       └── h_vol/
│           ├── indexes/
│           └── genomic.gff
├── workflow/
│   └── envs/
│       ├── fastp.yaml
│       ├── qc.yaml
│       └── rna_seq_processing.yaml
├── profile/
|   └── config.yaml
```

---

## Troubleshooting

### Common Issues
1. **Missing Conda Environments**:
   - Ensure you created the environments using the provided YAML files.
   - Verify with `conda env list`.

2. **Incorrect File Paths**:
   - Check that the paths in `samples.csv` and `config.yaml` are correct.

3. **Snakemake Errors**:
   - Use the `--printshellcmds` flag to debug:
     ```bash
     snakemake --use-conda --cores <number_of_cores> --printshellcmds
     ```

4. **Dependency Issues**:
   - If a tool is missing, recreate the environment:
     ```bash
     conda env create -f workflow/envs/<env_file>.yaml
     ```

