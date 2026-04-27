#!/bin/bash
#SBATCH --job-name=bismark2bedgraph
#SBATCH --mem=96G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=2-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

module load bowtie2/2.4.2
module load bismark/0.23.0

EXTRACT_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/methylation_extraction_post_mbias"
BEDGRAPH_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/bedgraph"

mkdir -p "$BEDGRAPH_DIR"

shopt -s nullglob

for sample_dir in "$EXTRACT_DIR"/*/ ; do
    sample=$(basename "$sample_dir")
    sample_out="$BEDGRAPH_DIR/$sample"

    mkdir -p "$sample_out"

    echo "Starting bismark2bedGraph for sample: $sample"
    echo "Input directory: $sample_dir"
    echo "Output directory: $sample_out"

    (
        cd "$sample_dir" || exit 1

        cpg_files=(CpG_*deduplicated*.txt CpG_*deduplicated*.txt.gz)
        chg_files=(CHG_*deduplicated*.txt CHG_*deduplicated*.txt.gz)
        chh_files=(CHH_*deduplicated*.txt CHH_*deduplicated*.txt.gz)

        if [ "${#cpg_files[@]}" -gt 0 ]; then
            echo "Processing CpG methylation calls for $sample"

            bismark2bedGraph \
                --ample_memory \
                --CX \
                --cutoff 5 \
                --output "$sample_out/bedgraph_output_cpg" \
                "${cpg_files[@]}"
        else
            echo "No CpG methylation extractor files found for $sample; skipping CpG."
        fi

        if [ "${#chg_files[@]}" -gt 0 ]; then
            echo "Processing CHG methylation calls for $sample"

            bismark2bedGraph \
                --ample_memory \
                --CX \
                --cutoff 5 \
                --output "$sample_out/bedgraph_output_chg" \
                "${chg_files[@]}"
        else
            echo "No CHG methylation extractor files found for $sample; skipping CHG."
        fi

        if [ "${#chh_files[@]}" -gt 0 ]; then
            echo "Processing CHH methylation calls for $sample"

            bismark2bedGraph \
                --ample_memory \
                --CX \
                --cutoff 5 \
                --output "$sample_out/bedgraph_output_chh" \
                "${chh_files[@]}"
        else
            echo "No CHH methylation extractor files found for $sample; skipping CHH."
        fi
    )

    echo "Creating uncompressed methylKit-compatible coverage files for sample: $sample"

    (
        cd "$sample_out" || exit 1

        cpg_cov=(*cpg*.bismark.cov.gz)
        chg_cov=(*chg*.bismark.cov.gz)
        chh_cov=(*chh*.bismark.cov.gz)

        if [ "${#cpg_cov[@]}" -gt 0 ]; then
            gunzip -c "${cpg_cov[0]}" > meth_CpG_cov_reads
        else
            echo "No CpG .bismark.cov.gz file found for $sample."
        fi

        if [ "${#chg_cov[@]}" -gt 0 ]; then
            gunzip -c "${chg_cov[0]}" > meth_CHG_cov_reads
        else
            echo "No CHG .bismark.cov.gz file found for $sample."
        fi

        if [ "${#chh_cov[@]}" -gt 0 ]; then
            gunzip -c "${chh_cov[0]}" > meth_CHH_cov_reads
        else
            echo "No CHH .bismark.cov.gz file found for $sample."
        fi
    )

    echo "Finished bismark2bedGraph for sample: $sample"
    echo
done

echo "All bismark2bedGraph processing complete."