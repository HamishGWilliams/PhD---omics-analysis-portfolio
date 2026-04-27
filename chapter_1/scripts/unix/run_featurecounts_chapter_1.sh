#!/bin/bash
#SBATCH --job-name=featureCounts_C1
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=12:00:00
#SBATCH --output=../../logs/outputs/featureCounts_C1_%j.out
#SBATCH --error=../../logs/errors/featureCounts_C1_%j.err

set -euo pipefail

# -----------------------------------------------------------------------------
# run_featurecounts_chapter_1.sh
#
# Runs featureCounts on Chapter 1 RNA-seq BAM files.
#
# Expected repository layout:
#   chapter_1/
#     data/
#       processed/              # input BAM files, or processed/bam/
#       external/                # annotation files, e.g. A_Equina.gtf
#     results/
#       tables/                  # featureCounts output table
#     logs/
#       outputs/                 # SLURM stdout logs
#       errors/                  # SLURM stderr logs
#     scripts/
#       unix/
#         run_featurecounts_chapter_1.sh
#
# Submit from the script directory or from the repository root using:
#   sbatch chapter_1/scripts/unix/run_featurecounts_chapter_1.sh
# -----------------------------------------------------------------------------

printf '\nRunning featureCounts for Chapter 1 RNA-seq BAM files\n'
printf 'Started at: %s\n' "$(date)"
printf 'SLURM job ID: %s\n\n' "${SLURM_JOB_ID:-not_running_under_slurm}"

# Load subread / featureCounts.
module load subread/2.0.2

# Resolve paths relative to this script, not relative to the submission directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHAPTER_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CHAPTER_DIR/.." && pwd)"

# Input locations.
# Prefer chapter_1/data/processed/bam if it exists; otherwise use chapter_1/data/processed.
if [[ -d "$CHAPTER_DIR/data/processed/bam" ]]; then
    BAM_DIR="$CHAPTER_DIR/data/processed/bam"
else
    BAM_DIR="$CHAPTER_DIR/data/processed"
fi

ANNOTATION="$CHAPTER_DIR/data/external/A_Equina.gtf"

# Output locations.
OUT_DIR="$CHAPTER_DIR/results/tables/featurecounts"
LOG_OUTPUT_DIR="$CHAPTER_DIR/logs/outputs"
LOG_ERROR_DIR="$CHAPTER_DIR/logs/errors"
OUT_FILE="$OUT_DIR/A_Equina_Counts.txt"

mkdir -p "$OUT_DIR" "$LOG_OUTPUT_DIR" "$LOG_ERROR_DIR"

printf 'Repository root: %s\n' "$REPO_ROOT"
printf 'Chapter directory: %s\n' "$CHAPTER_DIR"
printf 'BAM directory: %s\n' "$BAM_DIR"
printf 'Annotation file: %s\n' "$ANNOTATION"
printf 'Output file: %s\n\n' "$OUT_FILE"

# Check required inputs.
if [[ ! -d "$BAM_DIR" ]]; then
    printf 'ERROR: BAM directory does not exist: %s\n' "$BAM_DIR" >&2
    exit 1
fi

if [[ ! -f "$ANNOTATION" ]]; then
    printf 'ERROR: Annotation file not found: %s\n' "$ANNOTATION" >&2
    printf 'Place A_Equina.gtf in chapter_1/data/external/ or update ANNOTATION in this script.\n' >&2
    exit 1
fi

# Find BAM files directly inside BAM_DIR. Change -maxdepth if BAM files are nested.
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)

if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
    printf 'ERROR: No .bam files found in %s\n' "$BAM_DIR" >&2
    exit 1
fi

printf 'Found %s BAM files:\n' "${#BAM_FILES[@]}"
printf '  %s\n' "${BAM_FILES[@]}"
printf '\n'

# Run featureCounts.
# -p : paired-end reads
# -M : count multi-mapping reads
# -s 2 : reversely stranded library; change to 0 for unstranded or 1 for forward-stranded
# -T : number of threads; matched to SLURM -n 8
# -t exon : count exon features
# -g gene_id : summarise counts by gene_id attribute in the GTF
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

printf '\nfeatureCounts complete.\n'
printf 'Output written to: %s\n' "$OUT_FILE"
printf 'Summary written to: %s.summary\n' "$OUT_FILE"
printf 'Finished at: %s\n' "$(date)"
