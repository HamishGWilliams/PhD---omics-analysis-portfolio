#!/bin/bash
#SBATCH --job-name=name
#SBATCH --mem=48G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1/logs/errors/%x_%j.err

cd /uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1
