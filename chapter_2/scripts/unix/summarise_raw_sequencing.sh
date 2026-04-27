cd /uoa/scratch/users/r02hw22/Methylation_Analyses/Data
indir="./Raw"
per_sample_out="raw_sequencing_per_sample.tsv"
stats_out="raw_sequencing_summary_stats.tsv"

if [ ! -d "$indir" ]; then
  echo "Error: directory $indir not found" >&2
  exit 1
fi

# Per-sample table
printf "sample\tsample_id\treads_R1\treads_R2\tread_pairs\ttotal_reads\tmean_read_length_R1\tmean_read_length_R2\tmean_base_quality_R1\tmean_base_quality_R2\traw_bases_total\n" > "$per_sample_out"

for f in "$indir"/*; do
  [ -f "$f" ] || continue

  awk -v fallback_name="$(basename "$f")" '
    BEGIN {
      sample = ""
      in_row = 0
      n = 0
      row = 0
    }

    /<title>Download page for/ {
      if (match($0, /Sample_[^ ]+/)) {
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

        # Use the smaller of R1/R2 as the number of complete read pairs
        read_pairs = (r1_reads < r2_reads ? r1_reads : r2_reads)
        total_reads = r1_reads + r2_reads
        raw_bases_total = (r1_reads * meanlen[1]) + (r2_reads * meanlen[2])

        sample_id = sample
        sub(/^Sample_/, "", sample_id)
        sub(/_D$/, "", sample_id)

        printf "%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.1f\t%.1f\t%.1f\t%.1f\t%.0f\n",
          sample, sample_id,
          r1_reads, r2_reads, read_pairs, total_reads,
          meanlen[1], meanlen[2], meanq[1], meanq[2], raw_bases_total
      }
    }
  ' "$f" >> "$per_sample_out"
done

# Sort by sample_id, keeping header
{
  head -n 1 "$per_sample_out"
  tail -n +2 "$per_sample_out" | sort -V -k2,2
} > "${per_sample_out}.tmp"
mv "${per_sample_out}.tmp" "$per_sample_out"

# Helper to calculate summary stats for a numeric column
calc_stats() {
  local col="$1"
  local metric="$2"

  tail -n +2 "$per_sample_out" | cut -f"$col" | awk -v metric="$metric" '
    {
      x[NR] = $1
      sum += $1
    }
    END {
      if (NR == 0) exit

      n = NR

      # sort values
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
      mean = sum / n

      if (n % 2 == 1) {
        median = x[(n + 1) / 2]
      } else {
        median = (x[n / 2] + x[n / 2 + 1]) / 2
      }

      printf "%s\t%d\t%.0f\t%.2f\t%.2f\t%.0f\t%.0f\n",
        metric, n, sum, mean, median, min, max
    }
  '
}

# Summary statistics table
printf "metric\tn_samples\tsum\tmean\tmedian\tmin\tmax\n" > "$stats_out"
calc_stats 5  "read_pairs"      >> "$stats_out"
calc_stats 6  "total_reads"     >> "$stats_out"
calc_stats 11 "raw_bases_total" >> "$stats_out"

echo "Wrote:"
echo "  $per_sample_out"
echo "  $stats_out"
