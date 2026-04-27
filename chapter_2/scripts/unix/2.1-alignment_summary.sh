#!/bin/bash
#SBATCH --job-name=summarise_alignment
#SBATCH --mem=4G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=00:30:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

CHAPTER_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2"

BAM_DIR="$CHAPTER_DIR/data/processed/bam"
OUT_DIR="$CHAPTER_DIR/results/tables"

mkdir -p "$OUT_DIR"

OUT_FILE="$OUT_DIR/alignment_summary.csv"

echo "sample,sample_id,read_pairs_total,total_reads,unique_best_hit_pairs,mapping_efficiency_reported,mapping_efficiency_calculated" > "$OUT_FILE"

shopt -s nullglob

for sample_dir in "$BAM_DIR"/*/ ; do
    sample=$(basename "$sample_dir")

    # Convert e.g. Sample_3-3_D -> 3
    # Convert e.g. Sample_29-11_redo_D -> 29_11
    sample_id="$sample"
    sample_id="${sample_id#Sample_}"
    sample_id="${sample_id%_D}"
    sample_id="${sample_id/_redo/}"
    sample_id="${sample_id/-/_}"

    report_files=(
        "$sample_dir"/*_PE_report.txt
        "$sample_dir"/*_bismark_report.txt
    )

    if [ "${#report_files[@]}" -eq 0 ]; then
        echo "No Bismark alignment report found for sample: $sample" >&2
        continue
    fi

    if [ "${#report_files[@]}" -gt 1 ]; then
        echo "More than one report found for sample: $sample; using first:" >&2
        printf '%s\n' "${report_files[@]}" >&2
    fi

    report_file="${report_files[0]}"

    awk -v sample="$sample" -v sample_id="$sample_id" '
    BEGIN {
        OFS=","
        total = ""
        aligned = ""
        rate = ""
    }

    /^Sequence pairs analysed in total:/ {
        total = $NF
    }

    /^Number of paired-end alignments with a unique best hit:/ {
        aligned = $NF
    }

    /^Mapping efficiency:/ {
        rate = $NF
    }

    END {
        if (total != "" && aligned != "" && rate != "") {
            printf "%s,%s,%s,%s,%s,%s,%.2f%%\n", \
                sample, \
                sample_id, \
                total, \
                total * 2, \
                aligned, \
                rate, \
                (aligned / total) * 100
        } else {
            printf "WARNING: Could not parse expected fields from %s\n", FILENAME > "/dev/stderr"
        }
    }
    ' "$report_file" >> "$OUT_FILE"

done

echo "Alignment summary written to:"
echo "$OUT_FILE"