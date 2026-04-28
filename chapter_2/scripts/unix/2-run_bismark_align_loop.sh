#!/bin/bash
#SBATCH --job-name=bismark_align
#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

module load bowtie2/2.4.2
module load bismark/0.23.0

TRIMMED_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/trimmed"
GENOME_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/genome_index"
OUT_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/bam"

mkdir -p "$OUT_DIR"

cd "$TRIMMED_DIR"

for d in ./*/ ; do
    sample=$(basename "$d")
    mkdir -p "$OUT_DIR/$sample"

    (
        cd "$d" && \
        bismark \
            --multicore 8 \
            --genome "$GENOME_DIR" \
            --output_dir "$OUT_DIR/$sample" \
            -1 *_R1_001.fastq.gz \
            -2 *_R2_001.fastq.gz
    )
done

