#!/usr/bin/env Rscript
#SBATCH --job-name=regionalise_methylation
#SBATCH --mem=400G
#SBATCH --partition=uoa-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=21-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

# Load Required Libraries
library(methylKit) # for methylation analyses
library(plyranges) # for manipulating genomic ranges
library(GenomicRanges) # for converting and using files as GRanges
library(qqman) # for creating Q-Q and Manhattan plots
library(dplyr) # for data manipulation
library(ggplot2) # for plotting
library(data.table) # for data table manipulation
library(Rsamtools)  # For bgzip compression

# set working directory
setwd("/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results")

# Set the base directories
export_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results"

# Initialize log file
log_file <- file.path(export_dir, "process_log.txt")
cat("Process Log\n", file = log_file)

# Load the annotated genome
genome <- read_gff3("../../Data2/combined_annotations.gff3")

# Convert to DF
genome <- as.data.frame(genome)
# replace NA data with "NA"
genome[] <- lapply(genome, function(x) {
    # Check if the column is a factor and handle it accordingly
    if(is.factor(x)) {
        # Convert factor to character to avoid levels issues
        x <- as.character(x)
    }
    # Replace NA with "NA" (as character string)
    x[is.na(x)] <- "NA"
    return(x)
})
genome <- as(genome, "GRanges")

# Sample IDs and treatments for each experiment
exp1_sample_ids <- c("3", "9", "13", "15", "16", "17", "22", "27")
exp2_sample_ids <- c("6", "7", "8", "10", "14", "20", "21", "29_11", "32_18", "36_2")

exp1_treatment <- c(1, 0, 0, 1, 0, 0, 1, 1)
exp2_treatment <- c(0, 0, 1, 0, 1, 1, 1, 0, 0, 1)

contexts <- c("CpG", "CHG", "CHH")

# Create a list of experiments
experiments <- list(
  exp1 = list(
    sample_ids = exp1_sample_ids,
    treatments = exp1_treatment
  ),
  exp2 = list(
    sample_ids = exp2_sample_ids,
    treatments = exp2_treatment
  )
)

# Loop over each experiment and context
for (exp_name in names(experiments)) { # exp1 & exp2
  exp_data <- experiments[[exp_name]] # subsets experiments to selected exp
  sample_ids <- exp_data$sample_ids # extracts the sample ids
  treatments <- exp_data$treatments # extracts the treatment ids
  
  for (context in contexts) {
    # Construct the directory name
    dir_name <- paste0(exp_name, "_", context) # e.g. "exp1_CpG" which = directory name
    data_dir <- file.path(export_dir, dir_name)
	    # example: "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results/exp1_CpG"
    
    # Check if the directory exists (error handling)
    if (!dir.exists(data_dir)) {
      cat("Data directory", data_dir, "does not exist. Skipping.\n", file = log_file, append = TRUE)
      next
    }
    
    # Initialize a list to hold methylRawDB objects and treatments
    methyl_db_list <- list() # Empty DB list
    valid_treatments <- c() # Empty treatment vector
    
    # Loop over samples
    for (i in seq_along(sample_ids)) { # for each sample id (e.g., 3, 9 or 27 ...)
    # The function `seq_along(sample_ids)` generates a sequence of 
    # numbers from 1 to the length of the `sample_ids` vector.
      sample_id <- sample_ids[i] # Selects each sample id iteratively: e.g. 3
      treatment <- treatments[i] # Selects the appropriate treatment id too
      
      # Construct the file name
      file_name <- paste0(sample_id, ".txt.bgz") # e.g., "s.txt.bgz"
      file_path <- file.path(data_dir, file_name) # generates the appropriate file path
      # "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results/exp1_CpG/3.txt.bgz"
      
      # Check if the file exists (error checking)
      if (!file.exists(file_path)) {
        cat("File", file_path, "does not exist. Skipping sample", sample_id, "\n", file = log_file, append = TRUE)
        next
      }  
      
      # Read the methylation data
      methyl_db <- readMethylDB(file_path)
      # reads in the specific file into a methyl_db object
      
      # Add to the list
      methyl_db_list[[sample_id]] <- methyl_db
	      # adds the individual file to the list, which retains data as the 
	      # loop recurs
      
      # Add the treatment
      valid_treatments <- c(valid_treatments, treatment)
	}     # Iteratively adds the treatments to a treatmnets list.
    
    # Check if we have any valid methylation data
    if (length(methyl_db_list) == 0) {
      cat("No methylation data found for experiment", exp_name, "context", context, "\n", file = log_file, append = TRUE)
      next
    }
	    # Checks if the data has been successfully extracted
    
    # Create methylRawListDB
    methyl_raw_list_db <- methylRawListDB(
      methyl_db_list,
      treatment = valid_treatments
    )
	    # generates the methylRawDBList object with the treatments
    
    # Get the unique types in the 'genome' data
    types <- unique(genome$type)
    
    # Loop over each unique 'type'
    for (current_type in types) {
      # Subset 'genome' to regions of the current type
      genome_type <- genome[genome$type == current_type]
      
      # Check if there are any regions of this type
      if (length(genome_type) == 0) {
        cat("No regions found for type:", current_type, "\n", file = log_file, append = TRUE)
        next  # Skip to the next type
      }
      
      cat("Processing regions of type:", current_type, "\n", file = log_file, append = TRUE)
      
      # Loop over each sample
      for (i in seq_along(methyl_raw_list_db)) {
        # Get the methylRawDB object for the current sample
        methyl_raw_db <- methyl_raw_list_db[[i]] # selects each sample individually
        sample_id <- names(methyl_raw_list_db)[i] # captures sample id name
        
        # Log sample processing start
        cat("Processing sample:", sample_id, "for type:", current_type, "\n", file = log_file, append = TRUE)
        
        # Retrieve data using getData()
        sample_data <- getData(methyl_raw_db) 
	        # gets data in df format
        
        # Convert methylation data to GRanges
        methyl_gr <- GRanges(
          seqnames = sample_data$chr,
          ranges = IRanges(start = sample_data$start, end = sample_data$end),
          strand = sample_data$strand,
          coverage = sample_data$coverage,
          numCs = sample_data$numCs,
          numTs = sample_data$numTs
        )
	        # change df to GRanges format, extracting from the appropriate columns
        
        # Ensure chromosome names match, pruning sequences not in genome_type
        seqlevels(methyl_gr, pruning.mode = "coarse") <- seqlevels(genome_type)
		    ## ESSENTIAL PART!
        
        # Find overlaps between methylation sites and regions of current_type
        overlaps <- findOverlaps(
          query = genome_type,
          subject = methyl_gr,
          ignore.strand = TRUE  # Set to FALSE if strand-specific
        )
        
        # Log number of overlaps found
        num_overlaps <- length(overlaps)
        cat("Number of overlaps found for sample", sample_id, "and type", current_type, ":", num_overlaps, "\n", file = log_file, append = TRUE)
        
        # Save overlaps information to a file
        overlaps_info_file <- file.path(
          export_dir, 
          exp_name, context, current_type, 
          paste0(sample_id, "_", current_type, "_overlaps_info.tsv")
        )
	        # makes a new file: example - "3_Gene_overlaps_info.tsv"
        
        # Create directories if they don't exist
        dir.create(dirname(overlaps_info_file), recursive = TRUE, showWarnings = FALSE)
        
        overlaps_info <- data.frame(
          Region = as.character(seqnames(genome_type[queryHits(overlaps)])),
          Region_Start = start(genome_type[queryHits(overlaps)]),
          Region_End = end(genome_type[queryHits(overlaps)]),
          Methylation_Site = as.character(seqnames(methyl_gr[subjectHits(overlaps)])),
          Site_Position = start(methyl_gr[subjectHits(overlaps)]),
          Coverage = mcols(methyl_gr)$coverage[subjectHits(overlaps)]
        )
        
        write.table(
          overlaps_info,
          file = overlaps_info_file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
        
        cat("Overlaps information for sample", sample_id, "and type", current_type, "saved to", overlaps_info_file, "\n", file = log_file, append = TRUE)
        
        # Initialize data.table to store results for this sample and type
        results_dt <- data.table(
          chr = as.character(seqnames(genome_type)),
          start = start(genome_type),
          end = end(genome_type),
          strand = as.character(strand(genome_type)),
          coverage = 0,
          numCs = 0,
          numTs = 0
        )
        
        # If overlaps are found, aggregate counts
        if (num_overlaps > 0) {
          # Extract overlapping indices
          region_idx <- queryHits(overlaps)
          methyl_idx <- subjectHits(overlaps)
          
          # Create a data.table with methylation data and corresponding region IDs
          overlap_dt <- data.table(
            region_id = region_idx,
            coverage = mcols(methyl_gr)$coverage[methyl_idx],
            numCs = mcols(methyl_gr)$numCs[methyl_idx],
            numTs = mcols(methyl_gr)$numTs[methyl_idx]
          )
          
          # Aggregate counts over regions
          agg_dt <- overlap_dt[, .(
            coverage = sum(coverage),
            numCs = sum(numCs),
            numTs = sum(numTs)
          ), by = region_id]
          
          # Update results_dt with aggregated counts
          results_dt[agg_dt$region_id, coverage := agg_dt$coverage]
          results_dt[agg_dt$region_id, numCs := agg_dt$numCs]
          results_dt[agg_dt$region_id, numTs := agg_dt$numTs]
        } else {
          cat("No overlaps found for sample", sample_id, "and type", current_type, ". All counts set to zero.\n", file = log_file, append = TRUE)
        }
        
        # Save the Results to a File
        output_file <- file.path(
          export_dir, 
          exp_name, context, current_type, 
          paste0(sample_id, "_", current_type, "_aggregated_methylation_counts.txt")
        )
        
        # Create directories if they don't exist
        dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
        
        write.table(
          results_dt,
          file = output_file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
	        # Outputs aggregated regionalized data into its own directory per
	        # experiment, context and region type.
        
        cat("Aggregated methylation counts for sample", sample_id, "and type", current_type, "have been saved to", output_file, "\n", file = log_file, append = TRUE)
      }
    }
  }
}

cat("Processing completed.\n", file = log_file, append = TRUE)