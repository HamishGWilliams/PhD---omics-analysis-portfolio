#!/bin/bash
#SBATCH --job-name=copy_bam_files
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

SOURCE_DIR="/uoa/scratch/users/r02hw22/Methylation_Analyses/Data2/Trimmed"
DEST_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/bam"

mkdir -p "$DEST_DIR"

echo "Finding BAM files in: $SOURCE_DIR"
echo "Copying BAM files to: $DEST_DIR"
echo "Preserving directory structure."

find "$SOURCE_DIR" -type f -name "*.bam" -print0 | while IFS= read -r -d '' bam_file; do
    relative_path="${bam_file#$SOURCE_DIR/}"
    dest_path="$DEST_DIR/$relative_path"
    dest_subdir=$(dirname "$dest_path")

    mkdir -p "$dest_subdir"

    echo "Copying:"
    echo "  From: $bam_file"
    echo "  To:   $dest_path"

    cp -v "$bam_file" "$dest_path"
done

echo "BAM copy complete."