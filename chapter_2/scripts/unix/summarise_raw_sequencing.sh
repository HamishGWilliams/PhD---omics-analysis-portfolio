#!/bin/bash
#SBATCH --job-name=summarise_raw_sequencing
#SBATCH --mem=4G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=00:30:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

set -euo pipefail

CHAPTER_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2"

RAW_DIR="$CHAPTER_DIR/data/raw"
OUT_DIR="$CHAPTER_DIR/results/tables"

PER_SAMPLE_OUT="$OUT_DIR/raw_sequencing_per_sample.tsv"
STATS_OUT="$OUT_DIR/raw_sequencing_summary_stats.tsv"

mkdir -p "$OUT_DIR"

if [ ! -d "$RAW_DIR" ]; then
  echo "Error: raw sequencing directory not found:" >&2
  echo "$RAW_DIR" >&2
  exit 1
fi

# -------------------------------
# Per-sample table
# -------------------------------

printf "sample\tsample_id\treads_R1\treads_R2\tread_pairs\ttotal_reads\tmean_read_length_R1\tmean_read_length_R2\tmean_base_quality_R1\tmean_base_quality_R2\traw_bases_total\n" > "$PER_SAMPLE_OUT"

shopt -s nullglob

for f in "$RAW_DIR"/*; do
  [ -f "$f" ] || continue

  awk -v fallback_name="$(basename "$f")" '
    BEGIN {
      sample = ""
      in_row = 0
      n = 0
      row = 0
    }

    /<title>Download page for/ {
      if (match($0, /Sample_[^ <]+/)) {
        sample = substr($0, RSTART, RLENGTH)
      }
    }

    /<tr>/ {
      in_row = 1
      n = 0
      delete fields
      next
    }

    in_row && /<td>/ {
      line = $0
      gsub(/.*<td>/, "", line)
      gsub(/<\/td>.*/, "", line)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
      fields[++n] = line
      next
    }

    in_row && /<\/tr>/ {
      in_row = 0

      # Data rows have a FASTQ filename in column 1
      if (n > 0 && fields[1] ~ /\.fastq\.gz$/) {
        row++
        file[row]    = fields[1]
        reads[row]   = fields[2] + 0
        meanlen[row] = fields[4] + 0
        meanq[row]   = fields[8] + 0
      }
    }

    END {
      if (sample == "") {
        sample = fallback_name
      }

      if (row >= 2) {
        r1_reads = reads[1]
        r2_reads = reads[2]

        # Use the smaller R1/R2 count as the number of complete read pairs
        read_pairs = (r1_reads < r2_reads ? r1_reads : r2_reads)
        total_reads = r1_reads + r2_reads
        raw_bases_total = (r1_reads * meanlen[1]) + (r2_reads * meanlen[2])

        sample_id = sample
        sub(/^Sample_/, "", sample_id)
        sub(/_D$/, "", sample_id)
        gsub(/_redo/, "", sample_id)
        gsub(/-/, "_", sample_id)

        printf "%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.1f\t%.1f\t%.1f\t%.1f\t%.0f\n",
          sample, sample_id,
          r1_reads, r2_reads, read_pairs, total_reads,
          meanlen[1], meanlen[2], meanq[1], meanq[2], raw_bases_total
      } else {
        printf "WARNING: Could not parse two FASTQ rows from %s\n", FILENAME > "/dev/stderr"
      }
    }
  ' "$f" >> "$PER_SAMPLE_OUT"
done

# -------------------------------
# Sort by sample_id, keeping header
# -------------------------------

{
  head -n 1 "$PER_SAMPLE_OUT"
  tail -n +2 "$PER_SAMPLE_OUT" | sort -V -k2,2
} > "${PER_SAMPLE_OUT}.tmp"

mv "${PER_SAMPLE_OUT}.tmp" "$PER_SAMPLE_OUT"

# -------------------------------
# Helper function for summary stats
# -------------------------------

calc_stats() {
  local col="$1"
  local metric="$2"

  tail -n +2 "$PER_SAMPLE_OUT" | cut -f"$col" | awk -v metric="$metric" '
    {
      x[NR] = $1
      sum += $1
    }

    END {
      if (NR == 0) exit

      n = NR

      for (i = 1; i <= n; i++) {
        for (j = i + 1; j <= n; j++) {
          if (x[i] > x[j]) {
            tmp = x[i]
            x[i] = x[j]
            x[j] = tmp
          }
        }
      }

      min = x[1]
      max = x[n]
      range = max - min
      mean = sum / n

      if (n % 2 == 1) {
        median = x[(n + 1) / 2]
      } else {
        median = (x[n / 2] + x[n / 2 + 1]) / 2
      }

      printf "%s\t%d\t%.0f\t%.2f\t%.2f\t%.0f\t%.0f\t%.0f\n",
        metric, n, sum, mean, median, min, max, range
    }
  '
}

# -------------------------------
# Summary statistics table
# -------------------------------

printf "metric\tn_samples\tsum\tmean\tmedian\tmin\tmax\trange\n" > "$STATS_OUT"

calc_stats 5  "read_pairs"      >> "$STATS_OUT"
calc_stats 6  "total_reads"     >> "$STATS_OUT"
calc_stats 11 "raw_bases_total" >> "$STATS_OUT"

echo "Wrote:"
echo "  $PER_SAMPLE_OUT"
echo "  $STATS_OUT"