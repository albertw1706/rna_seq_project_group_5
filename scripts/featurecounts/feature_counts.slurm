#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="feature_counts"
#SBATCH --output=featureCounts.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --partition=pibu_el8

# Create the output directory
mkdir -p /data/users/awidjaja/rna_seq_project/feature_counts/new_feature_counts

# Set input directory
input_dir="/data/users/awidjaja/rna_seq_project/map/sorted_bam_and_index/"

# Set all the samples used
srr_accesions=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# Set output directory
output_dir="/data/users/awidjaja/rna_seq_project/feature_counts/new_feature_counts"

# Set annotation file directory
annotation_dir="/data/users/awidjaja/rna_seq_project/map"

# Empty list to store sample's BAM file names
list_bam=()

# Loop over to get all the sample's BAM file names
for srr in "${srr_accesions[@]}"; do
result="${input_dir}${srr}_sorted.bam "
list_bam+="$result"
done

# Run featurecounts
apptainer exec -B ${input_dir}:/input -B ${output_dir}:/output -B ${annotation_dir}:/annotations /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -p -T 4 -s 2 -Q 10 -a /annotations/Mus_musculus.GRCm39.113.gtf -o /output/summary_feature_counts.txt $list_bam