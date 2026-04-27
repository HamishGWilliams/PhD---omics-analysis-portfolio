awk '
BEGIN {
  OFS=","
  print "sample","sample_id","read_pairs_total","total_reads","unique_best_hit_pairs","mapping_efficiency_reported","mapping_efficiency_calculated"
}

/^Input files to be analysed/ {
  if (sample != "" && total != "" && aligned != "" && rate != "") {
    printf "%s,%s,%s,%s,%s,%s,%.2f%%\n", sample, sample_id, total, total*2, aligned, rate, (aligned/total)*100
  }

  split($0, a, "/")
  sample = a[length(a)]
  sub(/'\):$/, "", sample)

  sample_id = sample
  sub(/^Sample_/, "", sample_id)
  sub(/_D$/, "", sample_id)

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
  if (sample != "" && total != "" && aligned != "" && rate != "") {
    printf "%s,%s,%s,%s,%s,%s,%.2f%%\n", sample, sample_id, total, total*2, aligned, rate, (aligned/total)*100
  }
}
' 2-alignment-COMPLETE.out > alignment_summary.csv
