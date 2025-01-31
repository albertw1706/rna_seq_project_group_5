# Transcriptomic Analysis of Mouse Samples in Response to Toxoplasmosis

## Repository Overview
This GitHub repository is designated for the project that is a part of the RNA-seq course from the University of Bern. The project focuses on identifying differential gene expression in sequenced mouse samples comparing those infected with Toxoplasma gondii (toxoplasmosis) to control samples. The repository contains scripts and explanations of the workflow. 

## Dataset Description
The dataset utilized for this project was obtained from a subset of the dataset provided by Singhania et al. (2019). The data consisted of sequencing information in FASTQ file format with its SRA accessions from NCBI. The data was sequenced from blood and lung tissue samples of Mus musculus. The samples were derived from two groups, which are healthy control mice and mice infected with Toxoplasma gondii. The metadata of each sequencing data is presented in Table 1. 

| Sample    | Tissue  | Condition |
|-----------|---------|-----------|
| SRR7821918 | Lung    | Infected  |
| SRR7821919 | Lung    | Infected  |
| SRR7821920 | Lung    | Infected  |
| SRR7821921 | Lung    | Infected  |
| SRR7821922 | Lung    | Infected  |
| SRR7821937 | Lung    | Control   |
| SRR7821938 | Lung    | Control   |
| SRR7821939 | Lung    | Control   |
| SRR7821949 | Blood   | Infected  |
| SRR7821950 | Blood   | Infected  |
| SRR7821951 | Blood   | Infected  |
| SRR7821952 | Blood   | Infected  |
| SRR7821953 | Blood   | Infected  |
| SRR7821968 | Blood   | Control   |
| SRR7821969 | Blood   | Control   |
| SRR7821970 | Blood   | Control   |

Table 1. Sample Metadata Information

## Sequencing Data Quality Control and Mapping 
Start with running FastQC to check the quality of the sequencing data by running [qc.slurm](scripts/QC/qc.slurm).

Once all the reports from each samples were obtained, run [multi-qc.slurm](scripts/QC/multi-qc.slurm) to run MultiQC that makes all the html report files into one single html report file.

## Read Alignment and Gene Count Quantification

After checking the quality of the sequencing data, collect the reference genome and GTF file by running :
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz 
wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

Both files can also be obtained from the [Ensembl FTP website](https://www.ensembl.org/info/data/ftp/index.html) by downloading from the DNA FASTA section and Gene sets GTF section for the species Mus Musculus.



## RNA-Seq Exploratory Data Analysis

## Differential Gene Expression Analysis

## Overrepresentation Analysis

