#!/bin/bash
#SBATCH --partition uoa-compute
#SBATCH --mem-per-cpu 48G
#SBATCH -c 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

# sbatch /uoa/home/r02hw22/Equina_Methylation_Analysis/Scripts/3-run_loop_bismark_deduplicate.sh

module load  bowtie2/2.4.2
module load  bismark/0.23.0

cd /uoa/home/r02hw22/Equina_Methylation_Analysis/Data2/Trimmed

for d in ./*/ ; do
        (cd "$d" && deduplicate_bismark --bam *.bam );

                    done
