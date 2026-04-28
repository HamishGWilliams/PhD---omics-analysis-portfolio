#!/bin/bash
#SBATCH --job-name=bismark_extract_post_mbias
#SBATCH --mem=48G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=3-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

module load bowtie2/2.4.2
module load bismark/0.23.0

DEDUP_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/deduplicated_bam"
EXTRACT_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/methylation_extraction_post_mbias"

mkdir -p "$EXTRACT_DIR"

cd "$DEDUP_DIR" || exit 1

shopt -s nullglob

for d in ./*/ ; do
    sample=$(basename "$d")
    sample_out="$EXTRACT_DIR/$sample"

    mkdir -p "$sample_out"

    echo "Starting post-M-bias methylation extraction for sample: $sample"
    echo "Input directory: $d"
    echo "Output directory: $sample_out"

    (
        cd "$d" || exit 1

        bam_files=(*deduplicated.bam)

        if [ "${#bam_files[@]}" -eq 0 ]; then
            echo "No deduplicated BAM found for sample: $sample; skipping."
            exit 0
        fi

        if [ "${#bam_files[@]}" -gt 1 ]; then
            echo "More than one deduplicated BAM found for sample: $sample; skipping to avoid ambiguity."
            printf '%s\n' "${bam_files[@]}"
            exit 1
        fi

        bismark_methylation_extractor \
            -p \
            --no_overlap \
            --comprehensive \
            --report \
            --multicore 8 \
            --ignore 2 \
            --ignore_r2 3 \
            --ignore_3prime 3 \
            --ignore_3prime_r2 3 \
            --bedGraph \
            --output "$sample_out" \
            "${bam_files[0]}"
    )

    echo "Finished post-M-bias methylation extraction for sample: $sample"
    echo
done

echo "All post-M-bias methylation extraction jobs complete."