#!/bin/bash
#SBATCH --job-name=bismark_dedup
#SBATCH --partition=uoa-compute
#SBATCH --mem=48G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

module load bowtie2/2.4.2
module load bismark/0.23.0

BAM_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/bam"
DEDUP_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/deduplicated_bam"

mkdir -p "$DEDUP_DIR"

cd "$BAM_DIR"

for d in ./*/ ; do
    sample=$(basename "$d")
    sample_out="$DEDUP_DIR/$sample"

    mkdir -p "$sample_out"

    echo "Deduplicating BAM for sample: $sample"

    (
        cd "$d" || exit 1

        bam_file=(*_bismark_bt2_pe.bam)

        if [ ! -e "${bam_file[0]}" ]; then
            echo "No Bismark paired-end BAM found in $d; skipping."
            exit 0
        fi

        if [ "${#bam_file[@]}" -gt 1 ]; then
            echo "More than one Bismark BAM found in $d; skipping to avoid ambiguity."
            printf '%s\n' "${bam_file[@]}"
            exit 1
        fi

        deduplicate_bismark \
            --bam \
            --output_dir "$sample_out" \
            "${bam_file[0]}"
    )

    echo "Finished sample: $sample"
    echo
done