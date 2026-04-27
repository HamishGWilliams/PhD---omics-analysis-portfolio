#!/bin/bash
#SBATCH --job-name=move_out_files
#SBATCH --mem=8G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=02:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/logs/errors/%x_%j.err

SOURCE="/uoa/scratch/shared/sbs/Actinia omics/NEOF RNA seq"
DEST="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/logs/outputs"

echo "Finding .out files"
echo "Source: $SOURCE"
echo "Destination: $DEST"
echo "Started at: $(date)"

mkdir -p "$DEST"

find "$SOURCE" -type f -name "*.out" -print0 | while IFS= read -r -d '' out_file; do
    filename="$(basename "$out_file")"

    echo "Moving: $out_file -> $DEST/$filename"

    mv -n "$out_file" "$DEST/"
done

echo "Move complete"
echo "Finished at: $(date)"

echo "Number of .out files now in destination:"
find "$DEST" -maxdepth 1 -type f -name "*.out" | wc -l