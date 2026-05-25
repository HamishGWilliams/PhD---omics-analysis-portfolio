#!/bin/bash
#SBATCH --job-name=copy_PE_reports
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

SOURCE_DIR="/uoa/scratch/users/r02hw22/Methylation_Analyses/Data2/Trimmed"
DEST_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/bam"

mkdir -p "$DEST_DIR"

echo "Searching for *PE_report.txt files in sample directories under:"
echo "$SOURCE_DIR"
echo
echo "Copying reports into matching sample directories under:"
echo "$DEST_DIR"
echo

for sample_dir in "$SOURCE_DIR"/*/ ; do
    sample=$(basename "$sample_dir")
    dest_sample_dir="$DEST_DIR/$sample"

    mkdir -p "$dest_sample_dir"

    echo "Processing sample: $sample"

    find "$sample_dir" -maxdepth 1 -type f -name "*PE_report.txt" -print0 | while IFS= read -r -d '' report_file; do
        echo "Copying:"
        echo "  From: $report_file"
        echo "  To:   $dest_sample_dir/"
        cp -v "$report_file" "$dest_sample_dir/"
    done
done

echo
echo "PE report copy complete."