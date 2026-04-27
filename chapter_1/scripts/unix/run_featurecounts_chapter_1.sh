#!/bin/bash
#SBATCH --job-name=featureCounts_C1
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=12:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/logs/outputs/featureCounts_C1_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/logs/errors/featureCounts_C1_%j.err

set -euo pipefail

# Run this script from the chapter_1 directory:
#   cd /uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1
#   sbatch scripts/unix/run_featurecounts_chapter_1.sh

module load subread/2.0.2

CHAPTER_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1"
BAM_DIR="$CHAPTER_DIR/data/processed"
EXTRACTS_DIR="$CHAPTER_DIR/data/processed/extracts"
OUT_DIR="$CHAPTER_DIR/results/tables/featurecounts"
OUT_FILE="$OUT_DIR/A_Equina_Counts.txt"

mkdir -p "$OUT_DIR" "$CHAPTER_DIR/logs/outputs" "$CHAPTER_DIR/logs/errors"

# Choose the exon extract / annotation file.
# Edit this line if your file has a different name.
ANNOTATION="$EXTRACTS_DIR/A_Equina.gtf"

# If A_Equina.gtf is not present, use the first .gtf, .gff, .gff3, or .saf file in extracts/.
if [[ ! -f "$ANNOTATION" ]]; then
    ANNOTATION="$(find "$EXTRACTS_DIR" -maxdepth 1 -type f \( -name "*.gtf" -o -name "*.gff" -o -name "*.gff3" -o -name "*.saf" \) | sort | head -n 1)"
fi

if [[ -z "$ANNOTATION" || ! -f "$ANNOTATION" ]]; then
    echo "ERROR: No annotation/exon extract file found in: $EXTRACTS_DIR" >&2
    exit 1
fi

mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)

if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No BAM files found in: $BAM_DIR" >&2
    exit 1
fi

echo "Running featureCounts"
echo "BAM directory: $BAM_DIR"
echo "Annotation/extract file: $ANNOTATION"
echo "Output file: $OUT_FILE"
echo "Number of BAM files: ${#BAM_FILES[@]}"
echo "Started: $(date)"

if [[ "$ANNOTATION" == *.saf ]]; then
    featureCounts \
        -p \
        -M \
        -s 2 \
        -T 8 \
        -F SAF \
        -a "$ANNOTATION" \
        -o "$OUT_FILE" \
        "${BAM_FILES[@]}"
else
    featureCounts \
        -p \
        -M \
        -s 2 \
        -T 8 \
        -t exon \
        -g gene_id \
        -a "$ANNOTATION" \
        -o "$OUT_FILE" \
        "${BAM_FILES[@]}"
fi

echo "featureCounts complete"
echo "Finished: $(date)"
