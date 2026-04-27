#!/usr/bin/env Rscript
#SBATCH --job-name=methylation_analysis
#SBATCH --mem=200G
#SBATCH --partition=uoa-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

# Setup
# Load Packages
library(methylKit) # for methylation analyses  
library(plyranges) # for data manipulation with genomic ranges  
library(GenomicRanges) # for converting and using files as GRanges  
library(dplyr) # for data manipulation  
library(ggplot2) # for plotting  
library(data.table) # for data table manipulation

# initialise parameters for combinations of loop
experiments <- c("exp1","exp2")
contexts <- c("CpG", "CHG", "CHH")
regions <- c("dispersed_repeat",
             "downstream_region",
             "exon",
             "five_prime_UTR",
             "gene",
             "ncRNA_gene",
             "promoter",
             "three_prime_UTR"
             )

# Initialize an empty vector to keep track of created data frame names
created_df_names_exp1 <- c()
created_df_names_exp2 <- c()

combined_df_exp1 <- NULL
combined_df_exp2 <- NULL

# loop over each combination
for (experiment in experiments){
  for (context in contexts){
    for (region in regions){
      
      # build df name
      data_frame_name <- paste0(experiment, "_",
                                context,"_",
                                region)
     
      # build file name
      file_name <- paste0("Actinia_DMBs_d5_",
                          experiment, "_",
                          context, "_",
                          region,
                          "without_outlier_df_ann.txt")
      
      # Read the file into a temporary data frame
      df <- read.table(
        file = file_name,
        header = TRUE,          # set to FALSE if your data lacks header
        sep = "\t",             # change if you're using a different delimiter
        stringsAsFactors = FALSE
      )
      
      # Add the experiment, context & region columns
      df$experiment <- experiment
      df$context <- context
      df$region <- region
      
      # Clean data
      #df <- unique(df)
      
      # Assign this data frame to an object with the chosen name
      assign(data_frame_name, df)
      
      # rbind to a new df
      if (experiment == "exp1"){
        created_df_names_exp1 <- c(created_df_names_exp1, data_frame_name)
      } else if (experiment == "exp2"){
        created_df_names_exp2 <- c(created_df_names_exp2, data_frame_name)
      }
      
    }
  }
  # clean up the df at the end:
  rm(df)
}


# Now run a loop to perform unique on all files:
# Separate loop to apply unique() to each data frame
for (experiment in experiments) {
  for (context in contexts) {
    for (region in regions) {
      
      # Reconstruct the name of the data frame
      data_frame_name <- paste0(experiment, "_", context, "_", region)
      
      # Retrieve the data frame from your global environment
      df <- get(data_frame_name)
      
      # Apply unique()
      df_unique <- unique(df)
      
      # Optionally reassign it to replace the old one
      assign(data_frame_name, df_unique)
    }
  }
}

# AFTER the loop is done, combine them all:
combined_df_exp1 <- do.call(rbind, mget(created_df_names_exp1))
combined_df_exp2 <- do.call(rbind, mget(created_df_names_exp2))

# remove all the individual data frames to free up memory
rm(list = created_df_names_exp1)
rm(list = created_df_names_exp2)

# Experiment 1
# Change names of levels of 'region' to look nicer
combined_df_exp1$region <- as.factor(combined_df_exp1$region)
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "dispersed_repeat"] <- "Repeat"
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "downstream_region"] <- "Downstream"
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "exon"] <- "Exon"
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "five_prime_UTR"] <- "5` UTR"
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "gene"] <- "Gene"
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "ncRNA_gene"] <- "ncRNA"
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "promoter"] <- "Promoter"
levels(combined_df_exp1$region)[levels(combined_df_exp1$region) == "three_prime_UTR"] <- "3` UTR"

# Experiment 2
# Change names of levels of 'region' to look nicer
combined_df_exp2$region <- as.factor(combined_df_exp2$region)
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "dispersed_repeat"] <- "Repeat"
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "downstream_region"] <- "Downstream"
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "exon"] <- "Exon"
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "five_prime_UTR"] <- "5` UTR"
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "gene"] <- "Gene"
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "ncRNA_gene"] <- "ncRNA"
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "promoter"] <- "Promoter"
levels(combined_df_exp2$region)[levels(combined_df_exp2$region) == "three_prime_UTR"] <- "3` UTR"

# Save combined dataframes as .txt files
write.table(combined_df_exp1, "./exp1_DMRs.txt", sep = "\t")
write.table(combined_df_exp2, "./exp2_DMRs.txt", sep = "\t")