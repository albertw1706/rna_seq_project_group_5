#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="index"
#SBATCH --output=index.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH --partition=pibu_el8

module load apptainer

# Create and go to the TP directory (index will be used for mapping)
mkdir -p /data/users/awidjaja/rna_seq_project/map
cd /data/users/awidjaja/rna_seq_project/map

# Ftp download the reference genomes (unmasked) and annotations (for featureCounts)
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

# Unzip the reference genome and annotation files
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.113.gtf.gz

in_dir="/data/users/awidjaja/rna_seq_project/map"
out_dir="/data/users/awidjaja/rna_seq_project/map"

# Execute hisat2-build to index the reference genome
apptainer exec -B ${in_dir}:/input -B ${out_dir}:/output /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus_GRCm39_index