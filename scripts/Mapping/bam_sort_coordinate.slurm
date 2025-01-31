#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="QC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=30G
#SBATCH --partition=pibu_el8

SAMPLE=$1

#create and the results directory
mkdir -p /data/users/awidjaja/rna_seq_project/map/bam
mkdir -p /data/users/awidjaja/rna_seq_project/map/sorted_bam_and_index

# Set SAM to BAM conversion directories
in_convert_dir="/data/users/awidjaja/rna_seq_project/map/sam"
out_convert_dir="/data/users/awidjaja/rna_seq_project/map/bam"

# Execute samtools view to convert SAM to BAM
apptainer exec -B ${in_convert_dir}:/input -B ${out_convert_dir}:/output /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools view -S -b /input/${SAMPLE}.sam -o /output/${SAMPLE}.bam

# Set sorting directories
in_sort_dir="/data/users/awidjaja/rna_seq_project/map/bam"
out_sort_dir="/data/users/awidjaja/rna_seq_project/map/sorted_bam_and_index"

# Execute samtools sort to sort BAM files
apptainer exec -B ${in_sort_dir}:/input -B ${out_sort_dir}:/output /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools sort -o /output/${SAMPLE}_sorted.bam /input/${SAMPLE}.bam

# Set indexing directories
index_dir="/data/users/awidjaja/rna_seq_project/map/sorted_bam_and_index"

# Execute samtools index to index sorted BAM files
apptainer exec -B ${index_dir}:/input -B ${index_dir}:/output /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools index /input/${SAMPLE}_sorted.bam /output/${SAMPLE}_index.bai

