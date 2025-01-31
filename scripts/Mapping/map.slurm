#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="map"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --partition=pibu_el8

## Run this SLURM code with the bash file "run_map.sh" to loop through all samples and run in parallel

# Set sample name from the first argument in sbatch command
SAMPLE=$1

#create the results directory
mkdir -p cd /data/users/awidjaja/rna_seq_project/map/sam

# Set the input samples, input reference, and output SAM file directories
in_sample_dir="/data/users/awidjaja/rna_seq_project/input_seq"
in_reference_dir="/data/users/awidjaja/rna_seq_project/map"
out_dir="/data/users/awidjaja/rna_seq_project/map/sam"

# Make sure standard output for each samples are named correctly for quality check
exec >"${out_dir}/${SAMPLE}_summary_hisat2.out" 2>&1

# Execute hisat2-build to index the reference genome
apptainer exec \
    -B ${in_sample_dir}:/input_sample \
    -B ${in_reference_dir}:/input_reference \
    -B ${out_dir}:/output /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    hisat2 \
    -x /input_reference/Mus_musculus_GRCm39_index \
    -1 /input_sample/${SAMPLE}_1.fastq.gz \
    -2 /input_sample/${SAMPLE}_2.fastq.gz \
    -S /output/${SAMPLE}.sam \ 
    --rna-strandness RF


