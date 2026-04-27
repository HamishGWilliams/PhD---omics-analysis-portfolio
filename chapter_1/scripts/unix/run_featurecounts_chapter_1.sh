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
#       processed/
#         *.bam                  # input BAM files, or processed/bam/*.bam
#         extracts/              # exon annotation/extract file for featureCounts
#       processed/bam/           # optional preferred BAM input folder
#     results/
#       tables/featurecounts/    # featureCounts output table
#     logs/
#       outputs/                 # SLURM stdout logs
#       errors/                  # SLURM stderr logs
#     scripts/
#       unix/
#         run_featurecounts_chapter_1.sh
#
# Submit from the repository root using:
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

# Exon annotation / extract directory used by featureCounts.
EXTRACTS_DIR="$CHAPTER_DIR/data/processed/extracts"

# Optional: set ANNOTATION_BASENAME if the extracts directory contains more than
# one possible annotation file.
# Example:
#   ANNOTATION_BASENAME="A_Equina_exons.saf"
ANNOTATION_BASENAME=""

# Output locations.
OUT_DIR="$CHAPTER_DIR/results/tables/featurecounts"
LOG_OUTPUT_DIR="$CHAPTER_DIR/logs/outputs"
LOG_ERROR_DIR="$CHAPTER_DIR/logs/errors"
OUT_FILE="$OUT_DIR/A_Equina_Counts.txt"

mkdir -p "$OUT_DIR" "$LOG_OUTPUT_DIR" "$LOG_ERROR_DIR"

printf 'Repository root: %s\n' "$REPO_ROOT"
printf 'Chapter directory: %s\n' "$CHAPTER_DIR"
printf 'BAM directory: %s\n' "$BAM_DIR"
printf 'Extracts directory: %s\n' "$EXTRACTS_DIR"
printf 'Output file: %s\n\n' "$OUT_FILE"

# Check required inputs.
if [[ ! -d "$BAM_DIR" ]]; then
    printf 'ERROR: BAM directory does not exist: %s\n' "$BAM_DIR" >&2
    exit 1
fi

if [[ ! -d "$EXTRACTS_DIR" ]]; then
    printf 'ERROR: Exon extracts directory does not exist: %s\n' "$EXTRACTS_DIR" >&2
    exit 1
fi

# Find the featureCounts annotation/extract file.
# Supported:
#   .gtf / .gff / .gff3 : featureCounts default GTF/GFF mode
#   .saf                : featureCounts SAF mode using -F SAF
if [[ -n "$ANNOTATION_BASENAME" ]]; then
    ANNOTATION="$EXTRACTS_DIR/$ANNOTATION_BASENAME"
    if [[ ! -f "$ANNOTATION" ]]; then
        printf 'ERROR: ANNOTATION_BASENAME was set but file was not found: %s\n' "$ANNOTATION" >&2
        exit 1
    fi
else
    mapfile -t ANNOTATION_FILES < <(
        find "$EXTRACTS_DIR" -maxdepth 1 -type f \
            \( -name "*.gtf" -o -name "*.gff" -o -name "*.gff3" -o -name "*.saf" \) \
            | sort
    )

    if [[ "${#ANNOTATION_FILES[@]}" -eq 0 ]]; then
        printf 'ERROR: No .gtf, .gff, .gff3, or .saf annotation/extract file found in %s\n' "$EXTRACTS_DIR" >&2
        exit 1
    fi

    if [[ "${#ANNOTATION_FILES[@]}" -gt 1 ]]; then
        printf 'ERROR: More than one possible annotation/extract file found in %s\n' "$EXTRACTS_DIR" >&2
        printf 'Set ANNOTATION_BASENAME in this script to choose one. Candidates:\n' >&2
        printf '  %s\n' "${ANNOTATION_FILES[@]}" >&2
        exit 1
    fi

    ANNOTATION="${ANNOTATION_FILES[0]}"
fi

# Determine annotation format.
ANNOTATION_FORMAT_ARGS=()
case "$ANNOTATION" in
    *.saf)
        ANNOTATION_FORMAT_ARGS=(-F SAF)
        ;;
    *.gtf|*.gff|*.gff3)
        ANNOTATION_FORMAT_ARGS=()
        ;;
    *)
        printf 'ERROR: Unsupported annotation file extension: %s\n' "$ANNOTATION" >&2
        exit 1
        ;;
esac

printf 'Annotation/extract file: %s\n' "$ANNOTATION"
if [[ "${#ANNOTATION_FORMAT_ARGS[@]}" -gt 0 ]]; then
    printf 'Annotation format arguments: %s\n' "${ANNOTATION_FORMAT_ARGS[*]}"
else
    printf 'Annotation format arguments: default GTF/GFF mode\n'
fi
printf '\n'

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
# -t exon : count exon features for GTF/GFF input
# -g gene_id : summarise counts by gene_id attribute for GTF/GFF input
#
# For SAF input, featureCounts uses the GeneID column from the SAF file.
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

printf '\nfeatureCounts complete.\n'
printf 'Output written to: %s\n' "$OUT_FILE"
printf 'Summary written to: %s.summary\n' "$OUT_FILE"
printf 'Finished at: %s\n' "$(date)"
