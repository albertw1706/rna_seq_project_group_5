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

## Sequencing Data Quality Control and Genome Assembly
- Start with running FastQC to check the quality of the sequencing data by running [qc.slurm](scripts/QC/qc.slurm).

- Once all the reports from each samples were obtained, run [multi-qc.slurm](scripts/QC/multi-qc.slurm) to run MultiQC that makes all the html report files into one single html report file.

After checking the quality of the sequencing data, collect the reference genome and GTF file by running :

wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz 

wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

Both files can also be obtained from the [Ensembl FTP website](https://www.ensembl.org/info/data/ftp/index.html) by downloading from the DNA FASTA section and Gene sets GTF section for the species Mus Musculus.

- Once the reference genome is downloaded, index the reference genome by running [index.slurm](scripts/Mapping/index.slurm), then run genome assembly with [run_map.sh](scripts/Mapping/run_map.sh). 
(Note: [run_map.sh](scripts/Mapping/run_map.sh) will loop over all samples and run the [map.slurm](scripts/Mapping/map.slurm) script to make sure all samples ran in parallel when working in a cluster)

- Next, convert the SAM files to BAM format, sort and then index each BAM files to prepare the files for the next step by running the [run_bam_sort_coordinate.sh](scripts/Mapping/run_bam_sort_coordinate.sh). 
(Note: [run_bam_sort_coordinate.sh](scripts/Mapping/run_bam_sort_coordinate.sh) will loop over all samples and run the [bam_sort_coordinate.slurm](scripts/Mapping/bam_sort_coordinate.slurm) script to make sure all samples ran in parallel when working in a cluster)

## Read Alignment and Gene Count Quantification
The sorted and indexed BAM files produced from previous process were then used as input for identifying host genes present in the sample and quantify it by calculating reads that were aligned with the exons from the annotation file. This process can be done by running [featurecounts.slurm](scripts/QC/multi-qc.slurm). 

## Differential Gene Expression Analysis in R

- Once the text file from featurecounts is obtained, Differential Gene Expression Analysis can be done. Most of the process was done using DEseq2. The first process would be to do the RNA-Seq Exploratory Data Analysis and visualize it with a PCA plot to check for similar gene expression profiles or outliers. Next, identifying the genes that were upregulated or downregulated in response of the infection, which can be visualized with a volcano plot and heatmaps. Lastly, Overrepresentation Analysis was done by doing Gene Ontology (GO) enrichment analysis to identify overrepresented biological processes, molecular functions, and cellular components associated from the genes.
- All of the commands are available in [something.r](aa)


