#!/bin/bash
#SBATCH --job-name=repeatmodeler_repeatmasker
#SBATCH --mem=96G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=5-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

set -euo pipefail

# -------------------------------
# Modules
# -------------------------------

module load repeatmodeler/2.0.2
module load repeatmasker/4.1.4

# -------------------------------
# Project paths
# -------------------------------

CHAPTER_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2"

EXTERNAL_DIR="$CHAPTER_DIR/data/external"
REPEAT_DIR="$CHAPTER_DIR/data/processed/repeats"

LOG_DIR="$CHAPTER_DIR/logs"
mkdir -p "$LOG_DIR/outputs" "$LOG_DIR/errors"

mkdir -p "$REPEAT_DIR"

# -------------------------------
# Input files
# -------------------------------

GENOME_FASTA="$EXTERNAL_DIR/A_Equina.fa"
FLANK_GFF="$EXTERNAL_DIR/a_equina_with_1kb_flanks.gff3"

# -------------------------------
# Output files / directories
# -------------------------------

DB_DIR="$REPEAT_DIR/repeatmodeler_db"
RMODELER_DIR="$REPEAT_DIR/repeatmodeler_output"
RMASKER_DIR="$REPEAT_DIR/repeatmasker_output"

DB_NAME="$DB_DIR/actinia_equina_TE_db"

REPEAT_GFF_RAW="$REPEAT_DIR/A_Equina_repeats_raw.gff3"
REPEAT_GFF_STANDARDISED="$REPEAT_DIR/A_Equina_repeats_dispersed_repeat.gff3"

COMBINED_GFF="$EXTERNAL_DIR/combined_annotations.gff3"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# -------------------------------
# Checks
# -------------------------------

if [ ! -f "$GENOME_FASTA" ]; then
    echo "ERROR: Genome FASTA not found:"
    echo "$GENOME_FASTA"
    exit 1
fi

if [ ! -f "$FLANK_GFF" ]; then
    echo "ERROR: Flank-annotated GFF3 not found:"
    echo "$FLANK_GFF"
    echo "Run the flank annotation script first."
    exit 1
fi

mkdir -p "$DB_DIR" "$RMODELER_DIR" "$RMASKER_DIR"

echo "Genome FASTA: $GENOME_FASTA"
echo "Flank GFF3:    $FLANK_GFF"
echo "Repeat output: $REPEAT_DIR"
echo "Threads:       $THREADS"

# -------------------------------
# 1. Build RepeatModeler database
# -------------------------------

echo "Building RepeatModeler database..."

if [ ! -f "${DB_NAME}.nhr" ] && [ ! -f "${DB_NAME}.nin" ] && [ ! -f "${DB_NAME}.nsq" ]; then
    BuildDatabase \
        -name "$DB_NAME" \
        -engine ncbi \
        "$GENOME_FASTA"
else
    echo "RepeatModeler database already appears to exist; skipping BuildDatabase."
fi

# -------------------------------
# 2. Run RepeatModeler
# -------------------------------

echo "Running RepeatModeler..."

cd "$RMODELER_DIR"

if [ ! -f "$RMODELER_DIR/consensi.fa.classified" ]; then
    RepeatModeler \
        -database "$DB_NAME" \
        -pa "$THREADS"

    # RepeatModeler writes to an RM_* directory.
    # Copy the classified consensus library to a stable location.
    CONSENSI_FILE=$(find "$RMODELER_DIR" -type f -name "consensi.fa.classified" | head -n 1)

    if [ -z "$CONSENSI_FILE" ]; then
        echo "ERROR: RepeatModeler finished, but no consensi.fa.classified file was found."
        exit 1
    fi

    cp "$CONSENSI_FILE" "$RMODELER_DIR/consensi.fa.classified"
else
    echo "RepeatModeler consensus library already exists; skipping RepeatModeler."
fi

REPEAT_LIBRARY="$RMODELER_DIR/consensi.fa.classified"

if [ ! -f "$REPEAT_LIBRARY" ]; then
    echo "ERROR: Repeat library not found:"
    echo "$REPEAT_LIBRARY"
    exit 1
fi

# -------------------------------
# 3. Run RepeatMasker
# -------------------------------

echo "Running RepeatMasker..."

GENOME_BASENAME=$(basename "$GENOME_FASTA")

if [ ! -f "$RMASKER_DIR/${GENOME_BASENAME}.out" ]; then
    RepeatMasker \
        -pa "$THREADS" \
        -lib "$REPEAT_LIBRARY" \
        -gff \
        -dir "$RMASKER_DIR" \
        "$GENOME_FASTA"
else
    echo "RepeatMasker .out file already exists; skipping RepeatMasker."
fi

RM_OUT="$RMASKER_DIR/${GENOME_BASENAME}.out"

if [ ! -f "$RM_OUT" ]; then
    echo "ERROR: RepeatMasker output file not found:"
    echo "$RM_OUT"
    exit 1
fi

# -------------------------------
# 4. Convert RepeatMasker .out to GFF3
# -------------------------------

echo "Converting RepeatMasker .out to GFF3..."

rmOutToGFF3.pl "$RM_OUT" > "$REPEAT_GFF_RAW"

# -------------------------------
# 5. Standardise repeat feature labels
# -------------------------------
# This makes the repeat features easier to use downstream.
# The feature type in column 3 is set to "dispersed_repeat",
# which matches the region label used in the downstream methylation scripts.

echo "Standardising repeat features as dispersed_repeat..."

awk '
BEGIN { FS=OFS="\t"; repeat_n=0 }

{
    if ($0 ~ /^#/) {
        print
        next
    }

    if (NF < 9) {
        next
    }

    repeat_n++

    $2 = "RepeatMasker"
    $3 = "dispersed_repeat"

    if ($9 == "." || $9 == "") {
        $9 = "ID=repeat_" repeat_n
    } else if ($9 !~ /(^|;)ID=/) {
        $9 = "ID=repeat_" repeat_n ";" $9
    }

    print
}
' "$REPEAT_GFF_RAW" > "$REPEAT_GFF_STANDARDISED"

# -------------------------------
# 6. Merge flank-annotated genome GFF3 with repeat annotations
# -------------------------------

echo "Merging flank annotations and repeat annotations..."

cat "$FLANK_GFF" "$REPEAT_GFF_STANDARDISED" > "$COMBINED_GFF"

echo "Repeat annotation workflow complete."
echo
echo "RepeatModeler library:"
echo "$REPEAT_LIBRARY"
echo
echo "Raw repeat GFF3:"
echo "$REPEAT_GFF_RAW"
echo
echo "Standardised repeat GFF3:"
echo "$REPEAT_GFF_STANDARDISED"
echo
echo "Combined annotation GFF3:"
echo "$COMBINED_GFF"