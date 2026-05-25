#!/bin/bash
#SBATCH --mem 48G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=3-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

# sbatch /uoa/home/r02hw22/Equina_Methylation_Analysis/Scripts/run_bismark_extract_postmbias_loop.sh

module load  bowtie2/2.4.2
module load  bismark/0.23.0

cd /uoa/home/r02hw22/Equina_Methylation_Analysis/Data2/Trimmed/

for d in ./*/ ; do
        (cd "$d" && bismark_methylation_extractor -p --no_overlap --comprehensive --report --multicore 8 --ignore 2 --ignore_r2 3 --ignore_3prime 3 --ignore_3prime_r2 3 --bedGraph  *deduplicated.bam ) ;
                    done
