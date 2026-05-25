#!/bin/bash
#SBATCH --job-name=move_dedup_bams
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

BAM_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/bam"
DEDUP_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/deduplicated_bam"

mkdir -p "$DEDUP_DIR"

cd "$BAM_DIR"

for d in ./*/ ; do
    sample=$(basename "$d")
    sample_out="$DEDUP_DIR/$sample"

    mkdir -p "$sample_out"

    echo "Searching sample directory: $sample"

    find "$d" -type f -name "*_pe.deduplicated.bam" -print0 | while IFS= read -r -d '' bam_file; do
        echo "Moving:"
        echo "  From: $bam_file"
        echo "  To:   $sample_out/"

        mv -v "$bam_file" "$sample_out/"
    done
done

echo "Finished moving deduplicated BAM files."