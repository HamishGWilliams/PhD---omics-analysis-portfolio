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

CHAPTER_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1"
BAM_DIR="$CHAPTER_DIR/data/processed"
ANNOTATION_DIR="$CHAPTER_DIR/data/external"
OUT_DIR="$CHAPTER_DIR/results/tables/featurecounts"
OUT_FILE="$OUT_DIR/A_Equina_Counts.txt"

mkdir -p "$OUT_DIR" "$CHAPTER_DIR/logs/outputs" "$CHAPTER_DIR/logs/errors"

# Load featureCounts if the module system is available.
# If the module command is broken/unavailable, continue only if featureCounts is already on PATH.
if command -v module >/dev/null 2>&1; then
    module load subread/2.0.2 || true
fi

if ! command -v featureCounts >/dev/null 2>&1; then
    echo "ERROR: featureCounts is not available." >&2
    echo "Tried: module load subread/2.0.2" >&2
    echo "Check the module name on this cluster, or load subread before submitting the job." >&2
    exit 1
fi

# Choose the annotation file.
# Edit this line if your file has a different name.
ANNOTATION="$ANNOTATION_DIR/A_Equina.gtf"

# If A_Equina.gtf is not present, use the first .gtf, .gff, or .gff3 file in data/external/.
if [[ ! -f "$ANNOTATION" ]]; then
    ANNOTATION="$(find "$ANNOTATION_DIR" -maxdepth 1 -type f \
        \( -name "*.gtf" -o -name "*.gff" -o -name "*.gff3" \) \
        | sort | head -n 1)"
fi

if [[ -z "$ANNOTATION" || ! -f "$ANNOTATION" ]]; then
    echo "ERROR: No GTF/GFF annotation file found in: $ANNOTATION_DIR" >&2
    echo "Expected one of: *.gtf, *.gff, *.gff3" >&2
    echo "Files currently in data/external/:" >&2
    find "$ANNOTATION_DIR" -maxdepth 1 -type f -printf '  %f\n' 2>/dev/null | sort >&2 || true
    exit 1
fi

mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)

if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No BAM files found in: $BAM_DIR" >&2
    exit 1
fi

echo "Running featureCounts"
echo "BAM directory: $BAM_DIR"
echo "Annotation directory: $ANNOTATION_DIR"
echo "Annotation file: $ANNOTATION"
echo "Output file: $OUT_FILE"
echo "Number of BAM files: ${#BAM_FILES[@]}"
echo "Started: $(date)"

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

echo "featureCounts complete"
echo "Finished: $(date)"
