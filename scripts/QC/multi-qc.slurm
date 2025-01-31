#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="MultiQC"
#SBATCH --output=multiqc.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=1G
#SBATCH --partition=pibu_el8

# Load MultiQC tool
module load MultiQC/1.11

# Set input and output directory
in_dir="/data/users/awidjaja/rna_seq_project/MultiQC/input"
out_dir="/data/users/awidjaja/rna_seq_project/MultiQC/results"

# Create the directories needed (input and output) and go to the output directory
mkdir -p /data/users/awidjaja/rna_seq_project/MultiQC
mkdir -p $in_dir
mkdir -p $out_dir
cd /data/users/awidjaja/rna_seq_project/MultiQC/input

# List of the input samples
srr_accesions=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# For loop through the SRR accessions
for srr in "${srr_accesions[@]}"; do

# Copy all the FASTQC reports to the inout directory
cp /data/users/awidjaja/rna_seq_project/QC/${srr}/* $in_dir

done

# Run MultiQC
multiqc ${in_dir} -o ${out_dir}
