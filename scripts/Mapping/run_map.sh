# List of all sample's SRR accessions 
srr_accesions=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# Loop through the SRR accessions to run each SLURM script with sbatch
for srr in "${srr_accesions[@]}"; do
sbatch map.slurm $srr
done