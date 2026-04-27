#!/bin/bash
#SBATCH --job-name=copy_bam_files
#SBATCH --mem=96G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00

SOURCE="/uoa/scratch/shared/sbs/Actinia omics/NEOF RNA seq"
DEST="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/data/processed/bam"

echo "Starting BAM file copy"
echo "Source: $SOURCE"
echo "Destination: $DEST"
echo "Job ID: $SLURM_JOB_ID"
echo "Started at: $(date)"

mkdir -p "$DEST"

echo "Finding BAM files..."

find "$SOURCE" -type f -name "*.bam" -print0 | while IFS= read -r -d '' bam_file; do
    echo "Copying: $bam_file"
    rsync -avh --progress --partial --append-verify "$bam_file" "$DEST/"
done

echo "Copy complete"
echo "Finished at: $(date)"

echo "Copied BAM files in destination:"
find "$DEST" -maxdepth 1 -type f -name "*.bam" | wc -l

echo "Destination size:"
du -sh "$DEST"