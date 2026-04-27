#!/bin/bash
#SBATCH --mem 16G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

cd /uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2

module load  bowtie2/2.4.2
module load  bismark/0.23.0

bismark_genome_preparation --bowtie2 /uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/genome_index
