infile="2-alignment-COMPLETE.out"
per_sample_out="cytosine_context_counts_per_sample.tsv"
stats_out="cytosine_context_summary_stats.tsv"

if [ ! -f "$infile" ]; then
  echo "Error: $infile not found" >&2
  exit 1
fi

awk '
BEGIN {
  OFS="\t"
  print "sample","sample_id","meth_CpG","unmeth_CpG","total_CpG","meth_CHG","unmeth_CHG","total_CHG","meth_CHH","unmeth_CHH","total_CHH","total_C_analysed"
}

function print_row() {
  if (sample != "" &&
      mcpg != "" && ucpg != "" &&
      mchg != "" && uchg != "" &&
      mchh != "" && uchh != "") {

    tcpg = mcpg + ucpg
    tchg = mchg + uchg
    tchh = mchh + uchh

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
      sample, sample_id,
      mcpg, ucpg, tcpg,
      mchg, uchg, tchg,
      mchh, uchh, tchh,
      totalC
  }
}

/^Input files to be analysed/ {
  print_row()

  split($0, a, "/")
  sample = a[length(a)]
  gsub(/[^[:alnum:]_.-]+$/, "", sample)

  sample_id = sample
  sub(/^Sample_/, "", sample_id)
  sub(/_D$/, "", sample_id)

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
  print_row()
}
' "$infile" > "$per_sample_out"

{
  head -n 1 "$per_sample_out"
  tail -n +2 "$per_sample_out" | sort -V -k2,2
} > "${per_sample_out}.tmp"
mv "${per_sample_out}.tmp" "$per_sample_out"

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
  for (i = 1; i <= n; i++) sum += x[i]
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
}
END {
  if (n1 > 0) print_stats("CpG_total", cpg, n1)
  if (n2 > 0) print_stats("CHG_total", chg, n2)
  if (n3 > 0) print_stats("CHH_total", chh, n3)
}
' "$per_sample_out" > "$stats_out"

echo "Wrote:"
echo "  $per_sample_out"
echo "  $stats_out"