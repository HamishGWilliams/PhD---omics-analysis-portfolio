#!/bin/bash
#SBATCH --mem 6G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

#SBATCH --time=1-00:00:00

# sbatch /uoa/home/r02hw22/Equina_Methylation_Analysis/Scripts/run_plotm-bias3.sh

cd /uoa/home/r02hw22/Equina_Methylation_Analysis/Data2/Trimmed/Sample_3-3_D

awk '{f="file" NR; print $0 > f}' RS='================'  *.M-bias.txt
#creates 7 files, first is not needed. others need only lines 2-152

module load r/4.2.2

R
#install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)

#setwd("/uoa/home/s02df9/equina_meth_analysis/'NEOF Methyl seq'/Trimmed/Sample_3-3_D")
#due to spaces in the column names skip those and add later
r1_cpg = read.table("file2", header = F, skip = 2, nrows = 151)
r1_chg = read.table("file3", header = F, skip = 2, nrows = 151)
r1_chh = read.table("file4", header = F, skip = 2, nrows = 151)

r2_cpg = read.table("file5", header = F, skip = 2, nrows = 151)
r2_chg = read.table("file6", header = F, skip = 2, nrows = 151)
r2_chh = read.table("file7", header = F, skip = 2, nrows = 151)

all_r1 = rbind(r1_cpg, r1_chg, r1_chh)

all_r2 = rbind(r2_cpg, r2_chg, r2_chh)

colnames(all_r1) = c("position", "meth", "un", "pc_meth", "coverage")
colnames(all_r2) = c("position", "meth", "un", "pc_meth", "coverage")

all_r1$context = rep(c("cpg", "chg", "chh"), each = nrow(r1_cpg))
all_r1$read = "1"

all_r2$context = rep(c("cpg", "chg", "chh"), each = nrow(r2_cpg))
all_r2$read = "2"
all_data = rbind(all_r1, all_r2)


read_meth_pc = ggplot(all_data, aes(x = position,
                                 y = pc_meth,
                                group = context,
                                col = context)) +
  geom_line() +
  theme_classic() +
  ylab("Percentage methylated") +
  facet_wrap(~read, ncol = 1)

png(file = "pecent methylated by position per read 3 - pre_mbias.png")
read_meth_pc
dev.off()
quit(save = "no")


