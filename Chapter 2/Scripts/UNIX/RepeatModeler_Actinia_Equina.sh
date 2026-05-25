#!/bin/bash 
#SBATCH --mem 96G  
#SBATCH -N 1  
#SBATCH -n 8
#SBATCH --partition uoa-compute
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=h.williams.22@abdn.ac.uk  
#SBATCH --time=5-00:00:00

# sbatch /uoa/scratch/users/ro2hw22/Methylation_Analyses/Data2/RepeatModeler_Actinia_Equina.sh

# Load in modules
module load repeatmodeler/2.0.2
module load repeatmasker/4.1.4

# Ensure in correct Directory
cd /uoa/scratch/users/ro2hw22/Methylation_Analyses/Data2/

# NOTE: I ran each step individually in order to not repeat steps. Uncomment each step if you haven't performed it yet

# Build the DB for the repeatModeler
# BuildDatabase -name actinia_equina_TE_db -engine ncbi A_Equina.fa

# Run the RepeatModeler
# RepeatModeler -database actinia_equina_TE_db -pa 4

# annotate Genome with repeatMasker
RepeatMasker -pa 4 -lib ./RM_497924.WedOct300822272024/consensi.fa.classified -gff A_Equina.fa

# Use rmOutToGff3 to convert .fa.out file to .ggf3 format
rmOutToGFF3.pl A_Equina.fa.out > A_Equina_ith_TEs.gff3

# Merge with Existing annotations
cat genome_plus_upstream_and_downstream.gff3 A_Equina_with_TEs.gff3 > combined_annotations.gff3