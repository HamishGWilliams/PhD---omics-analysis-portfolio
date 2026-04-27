#!/bin/bash
#SBATCH --job-name=annotate_flanks
#SBATCH --mem=8G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=01:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

set -euo pipefail

GFF_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/external"

INPUT_GFF="$GFF_DIR/a_equina.gff3"
OUTPUT_GFF="$GFF_DIR/a_equina_with_1kb_flanks.gff3"

UPSTREAM=1000
DOWNSTREAM=1000

awk -v upstream="$UPSTREAM" -v downstream="$DOWNSTREAM" '
BEGIN { FS=OFS="\t" }

function get_attr(attrs, key,   n, i, f, pair) {
    n = split(attrs, f, ";")
    for (i = 1; i <= n; i++) {
        split(f[i], pair, "=")
        if (pair[1] == key) return pair[2]
    }
    return ""
}

{
    if ($0 ~ /^#/) {
        print
        next
    }

    # Print the original GFF3 feature first
    print

    if ($3 == "gene") {
        gene_id = get_attr($9, "ID")
        chr     = $1
        start   = $4
        end     = $5
        strand  = $7

        if (gene_id != "" && (strand == "+" || strand == "-")) {

            # -------------------------------
            # Promoter region
            # -------------------------------
            if (strand == "+") {
                p_start = start - upstream
                if (p_start < 1) p_start = 1
                p_end = start - 1
            } else if (strand == "-") {
                p_start = end + 1
                p_end = end + upstream
            }

            if (p_start <= p_end) {
                printf "%s\tcustom\tpromoter\t%d\t%d\t.\t%s\t.\tID=promoter_of_%s;Parent=%s\n", \
                    chr, p_start, p_end, strand, gene_id, gene_id
            }

            # -------------------------------
            # Downstream region
            # -------------------------------
            if (strand == "+") {
                d_start = end + 1
                d_end = end + downstream
            } else if (strand == "-") {
                d_start = start - downstream
                if (d_start < 1) d_start = 1
                d_end = start - 1
            }

            if (d_start <= d_end) {
                printf "%s\tcustom\tdownstream_region\t%d\t%d\t.\t%s\t.\tID=downstream_of_%s;Parent=%s\n", \
                    chr, d_start, d_end, strand, gene_id, gene_id
            }
        }
    }
}
' "$INPUT_GFF" > "$OUTPUT_GFF"

echo "Finished annotating 1 kb promoter and downstream regions."
echo "Input:  $INPUT_GFF"
echo "Output: $OUTPUT_GFF"