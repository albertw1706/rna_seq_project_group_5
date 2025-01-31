#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="QC"
#SBATCH --output=fastqc.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=1G
#SBATCH --partition=pibu_el8

# Set input directory
in_dir="/data/users/awidjaja/rna_seq_project/input_seq"

# Create and move to the input samples directory
mkdir -p $in_dir
cd $in_dir

# Link input files to the input directory
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/*

# List of all sample SRR accessions
srr_accesions=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# For loop through the SRR accessions
for srr in "${srr_accesions[@]}"; do

# Create and go to the output QC directory
mkdir -p /data/users/awidjaja/rna_seq_project/QC/${srr}

# Set the output QC directory
out_dir="/data/users/awidjaja/rna_seq_project/QC/${srr}"

# Run FASTQC in apptainer
apptainer exec \
    -B ${in_dir}:/input \
    -B ${out_dir}:/output /containers/apptainer/fastqc-0.12.1.sif fastqc \
    -t 1 \
    -o /output /input/${srr}_1.fastq.gz /input/${srr}_2.fastq.gz

done