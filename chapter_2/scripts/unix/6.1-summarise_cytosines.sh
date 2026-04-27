#!/bin/bash
#SBATCH --job-name=summarise_cytosines
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

EXTRACT_DIR="$CHAPTER_DIR/data/processed/methylation_extraction_post_mbias"
OUT_DIR="$CHAPTER_DIR/results/tables"

mkdir -p "$OUT_DIR"

PER_SAMPLE_OUT="$OUT_DIR/cytosine_context_counts_per_sample.tsv"
STATS_OUT="$OUT_DIR/cytosine_context_summary_stats.tsv"

echo -e "sample\tsample_id\tmeth_CpG\tunmeth_CpG\ttotal_CpG\tmeth_CHG\tunmeth_CHG\ttotal_CHG\tmeth_CHH\tunmeth_CHH\ttotal_CHH\ttotal_C_analysed" > "$PER_SAMPLE_OUT"

shopt -s nullglob

for sample_dir in "$EXTRACT_DIR"/*/ ; do
    sample=$(basename "$sample_dir")

    # Convert sample names to compact IDs:
    # Sample_3-3_D          -> 3
    # Sample_29-11_redo_D   -> 29_11
    # Sample_32-18_redo_D   -> 32_18
    sample_id="$sample"
    sample_id="${sample_id#Sample_}"
    sample_id="${sample_id%_D}"
    sample_id="${sample_id/_redo/}"
    sample_id="${sample_id/-/_}"

    report_files=(
        "$sample_dir"/*splitting_report.txt
        "$sample_dir"/*_splitting_report.txt
    )

    if [ "${#report_files[@]}" -eq 0 ]; then
        echo "No Bismark methylation extractor splitting report found for sample: $sample" >&2
        continue
    fi

    if [ "${#report_files[@]}" -gt 1 ]; then
        echo "More than one splitting report found for sample: $sample; using first:" >&2
        printf '%s\n' "${report_files[@]}" >&2
    fi

    report_file="${report_files[0]}"

    awk -v sample="$sample" -v sample_id="$sample_id" '
    BEGIN {
      OFS="\t"
      mcpg = ucpg = mchg = uchg = mchh = uchh = totalC = ""
    }

    /^Total number of C.s analysed:/ {
      totalC = $NF
    }

    /^Total methylated C.s in CpG context:/ {
      mcpg = $NF
    }

    /^Total methylated C.s in CHG context:/ {
      mchg = $NF
    }

    /^Total methylated C.s in CHH context:/ {
      mchh = $NF
    }

    /^Total unmethylated C.s in CpG context:/ {
      ucpg = $NF
    }

    /^Total unmethylated C.s in CHG context:/ {
      uchg = $NF
    }

    /^Total unmethylated C.s in CHH context:/ {
      uchh = $NF
    }

    END {
      if (mcpg != "" && ucpg != "" &&
          mchg != "" && uchg != "" &&
          mchh != "" && uchh != "") {

        tcpg = mcpg + ucpg
        tchg = mchg + uchg
        tchh = mchh + uchh

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
          sample, sample_id, \
          mcpg, ucpg, tcpg, \
          mchg, uchg, tchg, \
          mchh, uchh, tchh, \
          totalC
      } else {
        printf "WARNING: Could not parse expected cytosine fields from %s\n", FILENAME > "/dev/stderr"
      }
    }
    ' "$report_file" >> "$PER_SAMPLE_OUT"

done

# Sort by sample_id while preserving header
{
  head -n 1 "$PER_SAMPLE_OUT"
  tail -n +2 "$PER_SAMPLE_OUT" | sort -V -k2,2
} > "${PER_SAMPLE_OUT}.tmp"

mv "${PER_SAMPLE_OUT}.tmp" "$PER_SAMPLE_OUT"

# Generate summary statistics across samples
awk -F'\t' '
function sort_array(x, n,   i, j, tmp) {
  for (i = 1; i <= n; i++) {
    for (j = i + 1; j <= n; j++) {
      if (x[i] > x[j]) {
        tmp = x[i]
        x[i] = x[j]
        x[j] = tmp
      }
    }
  }
}

function print_stats(name, x, n,   sum, mean, median, min, max, range, i) {
  sum = 0

  for (i = 1; i <= n; i++) {
    sum += x[i]
  }

  sort_array(x, n)

  min = x[1]
  max = x[n]
  range = max - min
  mean = sum / n

  if (n % 2 == 1) {
    median = x[(n + 1) / 2]
  } else {
    median = (x[n / 2] + x[n / 2 + 1]) / 2
  }

  printf "%s\t%d\t%.0f\t%.0f\t%.0f\t%.2f\t%.2f\n", name, n, min, max, range, mean, median
}

BEGIN {
  OFS="\t"
  print "context","n_samples","min","max","range","mean","median"
}

NR > 1 {
  cpg[++n1] = $5
  chg[++n2] = $8
  chh[++n3] = $11
  total[++n4] = $12
}

END {
  if (n1 > 0) print_stats("CpG_total", cpg, n1)
  if (n2 > 0) print_stats("CHG_total", chg, n2)
  if (n3 > 0) print_stats("CHH_total", chh, n3)
  if (n4 > 0) print_stats("total_C_analysed", total, n4)
}
' "$PER_SAMPLE_OUT" > "$STATS_OUT"

echo "Wrote:"
echo "  $PER_SAMPLE_OUT"
echo "  $STATS_OUT"