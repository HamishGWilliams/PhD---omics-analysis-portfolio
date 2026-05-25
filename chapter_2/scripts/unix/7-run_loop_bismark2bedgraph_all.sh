#!/bin/bash
#SBATCH --job-name=bismark2bedgraph
#SBATCH --mem=96G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
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
    sample_dir="${sample_dir%/}"
    sample=$(basename "$sample_dir")
    sample_out="$BEDGRAPH_DIR/$sample"

    mkdir -p "$sample_out"

    echo "============================================================"
    echo "Starting bismark2bedGraph for sample: $sample"
    echo "Input directory: $sample_dir"
    echo "Output directory: $sample_out"
    echo "============================================================"

    cpg_files=("$sample_dir"/CpG_*deduplicated*.txt)
    chg_files=("$sample_dir"/CHG_*deduplicated*.txt)
    chh_files=("$sample_dir"/CHH_*deduplicated*.txt)

    if [ "${#cpg_files[@]}" -eq 0 ]; then
        cpg_files=("$sample_dir"/CpG_*deduplicated*.txt.gz)
    fi

    if [ "${#chg_files[@]}" -eq 0 ]; then
        chg_files=("$sample_dir"/CHG_*deduplicated*.txt.gz)
    fi

    if [ "${#chh_files[@]}" -eq 0 ]; then
        chh_files=("$sample_dir"/CHH_*deduplicated*.txt.gz)
    fi

    (
        cd "$sample_out" || exit 1

        echo "Working from output directory:"
        pwd

        if [ "${#cpg_files[@]}" -gt 0 ]; then
            echo "Processing CpG methylation calls for $sample"

            cpg_links=()

            for file in "${cpg_files[@]}"; do
                link_name=$(basename "$file")
                ln -sf "$file" "$link_name"
                cpg_links+=("$link_name")
            done

            printf 'CpG input symlinks:\n'
            printf '  %s\n' "${cpg_links[@]}"

            rm -f bedgraph_output_cpg.gz
            rm -f bedgraph_output_cpg.gz.bismark.cov.gz

            bismark2bedGraph \
                --scaffolds \
                --CX \
                --cutoff 5 \
                --buffer_size 20G \
                --output "bedgraph_output_cpg" \
                "${cpg_links[@]}"
        else
            echo "No CpG methylation extractor files found for $sample; skipping CpG."
        fi

        if [ "${#chg_files[@]}" -gt 0 ]; then
            echo "Processing CHG methylation calls for $sample"

            chg_links=()

            for file in "${chg_files[@]}"; do
                link_name=$(basename "$file")
                ln -sf "$file" "$link_name"
                chg_links+=("$link_name")
            done

            printf 'CHG input symlinks:\n'
            printf '  %s\n' "${chg_links[@]}"

            rm -f bedgraph_output_chg.gz
            rm -f bedgraph_output_chg.gz.bismark.cov.gz

            bismark2bedGraph \
                --scaffolds \
                --CX \
                --cutoff 5 \
                --buffer_size 20G \
                --output "bedgraph_output_chg" \
                "${chg_links[@]}"
        else
            echo "No CHG methylation extractor files found for $sample; skipping CHG."
        fi

        if [ "${#chh_files[@]}" -gt 0 ]; then
            echo "Processing CHH methylation calls for $sample"

            chh_links=()

            for file in "${chh_files[@]}"; do
                link_name=$(basename "$file")
                ln -sf "$file" "$link_name"
                chh_links+=("$link_name")
            done

            printf 'CHH input symlinks:\n'
            printf '  %s\n' "${chh_links[@]}"

            rm -f bedgraph_output_chh.gz
            rm -f bedgraph_output_chh.gz.bismark.cov.gz

            bismark2bedGraph \
                --scaffolds \
                --CX \
                --cutoff 5 \
                --buffer_size 20G \
                --output "bedgraph_output_chh" \
                "${chh_links[@]}"
        else
            echo "No CHH methylation extractor files found for $sample; skipping CHH."
        fi

        echo "Creating uncompressed methylKit-compatible coverage files for sample: $sample"

        cpg_cov=(bedgraph_output_cpg*.bismark.cov.gz)
        chg_cov=(bedgraph_output_chg*.bismark.cov.gz)
        chh_cov=(bedgraph_output_chh*.bismark.cov.gz)

        if [ "${#cpg_cov[@]}" -gt 0 ]; then
            echo "Creating meth_CpG_cov_reads from ${cpg_cov[0]}"
            gunzip -c "${cpg_cov[0]}" > meth_CpG_cov_reads
        else
            echo "No CpG .bismark.cov.gz file found for $sample."
        fi

        if [ "${#chg_cov[@]}" -gt 0 ]; then
            echo "Creating meth_CHG_cov_reads from ${chg_cov[0]}"
            gunzip -c "${chg_cov[0]}" > meth_CHG_cov_reads
        else
            echo "No CHG .bismark.cov.gz file found for $sample."
        fi

        if [ "${#chh_cov[@]}" -gt 0 ]; then
            echo "Creating meth_CHH_cov_reads from ${chh_cov[0]}"
            gunzip -c "${chh_cov[0]}" > meth_CHH_cov_reads
        else
            echo "No CHH .bismark.cov.gz file found for $sample."
        fi

        echo "Removing temporary input symlinks for $sample"

        find . -maxdepth 1 -type l -name "CpG_*deduplicated*.txt" -delete
        find . -maxdepth 1 -type l -name "CHG_*deduplicated*.txt" -delete
        find . -maxdepth 1 -type l -name "CHH_*deduplicated*.txt" -delete
        find . -maxdepth 1 -type l -name "CpG_*deduplicated*.txt.gz" -delete
        find . -maxdepth 1 -type l -name "CHG_*deduplicated*.txt.gz" -delete
        find . -maxdepth 1 -type l -name "CHH_*deduplicated*.txt.gz" -delete
    )

    echo "Finished bismark2bedGraph for sample: $sample"
    echo
done

echo "All bismark2bedGraph processing complete."