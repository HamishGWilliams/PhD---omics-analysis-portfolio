#!/bin/bash
#SBATCH --job-name=copy_trimmed_emseq
#SBATCH --mem=96G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00

SOURCE="/uoa/scratch/users/r02hw22/Methylation_Analyses/Data2/Trimmed/"
DEST="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/trimmed/"

echo "Starting RNA-seq trimmed file copy"
echo "Source: $SOURCE"
echo "Destination: $DEST"
echo "Job ID: $SLURM_JOB_ID"
echo "Started at: $(date)"

mkdir -p "$DEST"

rsync -avh --progress --partial --append-verify \
  "$SOURCE" \
  "$DEST"

echo "Copy complete"
echo "Finished at: $(date)"

echo "Source size:"
du -sh "$SOURCE"

echo "Destination size:"
du -sh "$DEST"

echo "Source file count:"
find "$SOURCE" -type f | wc -l

echo "Destination file count:"
find "$DEST" -type f | wc -l