#!/bin/bash
#SBATCH --job-name=move_bedgraph
#SBATCH --mem=96G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

SOURCE_BASE="/uoa/scratch/users/r02hw22/Methylation_Analyses/Data2/Trimmed"
DEST_BASE="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/bedgraph"

for sample_dir in "$SOURCE_BASE"/Sample_*; do
    [ -d "$sample_dir" ] || continue

    sample_name="$(basename "$sample_dir")"
    dest_dir="$DEST_BASE/$sample_name"

    echo "Processing $sample_name"

    # Create matching target directory if it does not already exist
    mkdir -p "$dest_dir"

    find "$sample_dir" -maxdepth 1 -type f \( \
        -name "bedgraph_output_*" -o \
        -name "meth_*" \
    \) -exec cp -v {} "$dest_dir"/ \;
done

echo "Copy complete."