#!/usr/bin/env Rscript
#SBATCH --job-name=methylation_analysis
#SBATCH --mem=200G
#SBATCH --partition=uoa-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=21-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

# Setup
# Load Packages
library(methylKit) # for methylation analyses  
library(plyranges) # for data manipulation with genomic ranges  
library(GenomicRanges) # for converting and using files as GRanges  
library(qqman) # for Q-Q and Manhattan plots  
library(dplyr) # for data manipulation  
library(ggplot2) # for plotting  
library(data.table) # for data table manipulation
# library(ggforce) # For drawing ellipses
# library(concaveman) # for adding concavity to hull

# Set base_path
base_path <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results"

# change to base_path
setwd(base_path)

# Generate objects to operate loop:
experiments <- c("exp1", "exp2")
# contexts <- c("CpG","CHG") # removing CHH because the jobs run out of memory

	# Running just for CHH with 200G memory to try and brute force the analyses (28th April 2025)

contexts <- c("CHH")

# Load in additional data files for analyses:
genome <- read_gff3("../../Data2/combined_annotations.gff3")
sample_metadata <- read.table("./sample_metadata.txt", header = T)
sample_metadata$SAMPLE_ID <- as.factor(sample_metadata$SAMPLE_ID)

# TESTING OPTIONS: COMMENT OUT WHEN TESTING IS COMPLETE!
# experiment <- "exp2"
# context <- "CHH"

# Commencing loop execution ----
for (experiment in experiments){ # start loop for each experiment
  for (context in contexts){ # start looping through each context too...
    if (experiment == "exp1") { # used to help generate correct list of files and data
      
      # Ensure objects used are empty
      file.list <- NULL
      sample.id <- NULL
      treatment <- NULL
      
      # build correct file list
      file.list = list( # paste0("../../Data2/Trimmed/Sample_3-3_D/meth_", context, "_cov_reads"), #  outlier
                       paste0("../../Data2/Trimmed/Sample_9-9_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_13-13_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_15-15_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_16-16_D/meth_", context, "_cov_reads"),
                        # paste0("../../Data2/Trimmed/Sample_17-17_D/meth_", context, "_cov_reads"), # outlier
                       paste0("../../Data2/Trimmed/Sample_22-22_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_27-27_D/meth_", context, "_cov_reads")
      )

      # build the sample.id
      sample.id = list("9","13","15","16","22","27") # remove the 1st and 6th values when removing outlier
      
      # build the treatment
      treatment = c(0, 0, 1, 0, 1, 1) # remove the 1st and 6th values when removing outlier
      
      
    } else if (experiment == "exp2") { # kind of redundant for exp1 script since exp2 isnt included in loop
      
      # build correct file list
      file.list = list(paste0("../../Data2/Trimmed/Sample_6-6_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_7-7_D/meth_", context, "_cov_reads"),
                        # paste0("../../Data2/Trimmed/Sample_8-8_D/meth_", context, "_cov_reads"), # outlier
                       paste0("../../Data2/Trimmed/Sample_10-10_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_14-14_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_20-20_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_21-21_D/meth_", context, "_cov_reads"),
                       paste0("../../Data2/Trimmed/Sample_29-11_redo_D/meth_", context, "_cov_reads"),
                        # paste0("../../Data2/Trimmed/Sample_32-18_redo_D/meth_", context, "_cov_reads"), ' outlier
                       paste0("../../Data2/Trimmed/Sample_36-2_redo_D/meth_", context, "_cov_reads")
      )
      
      # build exp2 sample.id list
      sample.id = list("6", "7", "10", "14", "20", "21", "11", "2") # removed the 3rd (8) and 9th (18) samples
      
      # build exp2 treatment list, in order of samples
      treatment = c(0, 0, 0, 1, 1, 1, 0, 1) # removed the 3rd (1) AND 9TH (0) samples
    
      
    } # End of the if & else if loop for making the correct experiment data objects
    
    # Back to the second layer of the loop (Experiments > Contexts) 
    
    # Read the methylation data ----
    act_data <- NULL
    tryCatch({ 
      act_data <- methRead(
        file.list,
        sample.id = sample.id,
        assembly = "actinia",
        pipeline = "bismarkCoverage",
        treatment = treatment,
        context = context,
        dbtype = "tabix",
        mincov = 0
      )
    }, error = function(e) {
      message("Error reading methylation data for experiment: ", experiment, ", context: ", context)
      message("Error message: ", e$message)
    })
    
    # If methRead failed, no point in continuing to next steps
    if (is.null(act_data)) {
      # Clean up and continue to next iteration
      rm(list = setdiff(ls(), c("experiments", "contexts", "experiment", "context", "base_path","genome", "sample_metadata")))
      gc()
      next
    }
    
    # Plot methylation statistics ----
    tryCatch({
      filename <- paste0("raw_methylation_plots_act_", context, "_", experiment, "_removing_outliers",".png")
      png(filename, width = 12, height = 8, units = "in", res = 300)
      lapply(act_data, getMethylationStats, plot = TRUE)
      dev.off()
    }, error = function(e) {
      message("Error plotting methylation stats for ", experiment, ", ", context, ": ", e$message)
      dev.off()
    })
    
    # Plot raw coverage statistics ----
    tryCatch({
      filename <- paste0("raw_coverage_plots_act_", context, "_",  experiment, "_removing_outliers",".png")
      png(filename, width = 20, height = 10, units = "in", res = 300)
      lapply(act_data, getCoverageStats, plot = TRUE)
      dev.off()
    }, error = function(e) {
      message("Error plotting coverage stats for ", experiment, ", ", context, ": ", e$message)
      dev.off()
    })
    
    # Filter for coverage ----
    act_data_f <- NULL
    tryCatch({
      act_data_f <- filterByCoverage(
        act_data,
        hi.perc = 99.9,
        lo.count = 5,
        lo.perc = NULL,
        suffix = paste0(context, "_f")
      )
    }, error = function(e) {
      message("Error filtering coverage for ", experiment, ", ", context, ": ", e$message)
    })
    
    if (is.null(act_data_f)) {
      rm(list = setdiff(ls(), c("experiments", "contexts", "experiment", "context", "base_path","genome","sample_metadata")))
      gc()
      next
    }
    
    # Normalize coverage ----
    act_data_f_norm <- NULL
    tryCatch({
      act_data_f_norm <- normalizeCoverage(act_data_f, method = "median")
    }, error = function(e) {
      message("Error normalizing coverage for ", experiment, ", ", context,": ", e$message)
    })
    
    if (is.null(act_data_f_norm)) {
      rm(list = setdiff(ls(), c("experiments", "contexts", "experiment", "context", "base_path","genome","sample_metadata")))
      gc()
      next
    }
    
    # Set data.table option
    options(datatable.allow.cartesian = F) # Ensure cartesians are no generated in uniting samples
    
    # Unite methylation results ----
    act_data_normed_fu <- NULL
    tryCatch({
      act_data_normed_fu <- methylKit::unite(
        act_data_f_norm,
        destrand = FALSE,
        min.per.group = 3L,
        suffix = paste0(context, "_normed_fu3")
      )
    }, error = function(e) {
      message("Error uniting methylation data for ", experiment, ", ", context, ": ", e$message)
    })
    
    if (is.null(act_data_normed_fu)) {
      rm(list = setdiff(ls(), c("experiments", "contexts", "experiment", "context", "base_path","genome","sample_metadata")))
      gc()
      next
    }	
    
    ## Generate Clustering Plots ----
    tryCatch({
    filename <- paste0("clustering_plot", context, "_", experiment, "_removing_outliers",".png")
    png(filename, width = 18, height = 12, units = "cm", res = 300)
    clusterSamples(act_data_normed_fu,
                   filterByQuantile = F,
                   sd.threshold = 0.1,
                   dist="correlation",
                   method="ward.D",
                   plot=T)
    dev.off()
    }, error = function(e) {
      message("Error normalizing coverage for ", experiment, ", ", context,": ", e$message)
    })
    
    ## PCA plots ----
    if (experiment == "exp1"){
    tryCatch({
      # Perform PCA Analysis
      pca_results <- PCASamples(act_data_normed_fu, comp = c(1, 2), obj.return = T)
      pca_data <- as.data.frame(pca_results$x)
      pca_data$Sample <- rownames(pca_data)          
      pca_data$Sample <- as.factor(pca_data$Sample)
      
      # Calculate variance explained
      sdev <- pca_results$sdev
      variance <- sdev^2
      total_variance <- sum(variance)
      explained_variance <- variance / total_variance
      
      # Join with metadata
      pca_data <- pca_data %>%
        left_join(sample_metadata, by = c("Sample" = "SAMPLE_ID"))
      
      # Separate groups
      data_group_control <- pca_data %>% filter(Group == "X") # change this for exp2
      data_group_treatment <- pca_data %>% filter(Group == "A") # change this for exp2
      
      # Inner join to create *_control and *_treatment columns
      arrow_data <- inner_join(data_group_control, data_group_treatment, by = "Clone", suffix = c("_control", "_treatment")) %>%
        mutate(
          mid_x = (PC1_control + PC1_treatment) / 2,
          mid_y = (PC2_control + PC2_treatment) / 2
        )
      
      # Generate Hull data:
      df_hull <- pca_data %>%
        group_by(Group) %>%
        slice(chull(PC1, PC2))  # slice() keeps only hull points for each cluster
      
      
      combined_pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point(aes(color = Group), size = 5) +
        # geom_mark_hull(concavity = 5, expand = 0, radius = 0, aes(fill=Group, colour = Group)) +
       #  stat_ellipse(aes(fill = Group), 
       #               type = "t",      # t-distribution based
       #               alpha = 0.2,     # semi-transparent fill
       #               geom = "polygon", 
       #               color = NA) +    # no outline, color is guided by "fill" above
        
        geom_polygon(
          data = df_hull, 
          aes(x = PC1, y = PC2, fill = Group), 
          alpha = 0.2, 
          color = NA
        ) +
        
        geom_segment(
          data = arrow_data,
          aes(x = PC1_control, y = PC2_control, xend = PC1_treatment, yend = PC2_treatment),
          arrow = arrow(length = unit(0.2, "inches")),
          color = "black"
        ) +
        geom_label(data = arrow_data, aes(x = mid_x, y = mid_y, label = Clone), size = 3) +
        scale_color_manual(values = c("A" = "red", "X" = "blue")) +
        scale_fill_manual(values = c("A" = "red", "X" = "blue")) +
        theme_bw() +
        labs(
          title = "PCA Plot: PCs 1-2",
          x = paste0("PC", 1, " (", round(explained_variance[1] * 100, 1), "%)"),
          y = paste0("PC", 2, " (", round(explained_variance[2] * 100, 1), "%)"),
          color = "Group"
        )
      
      filename <- paste0("PCA_", 1, "-", 2, "_plot_", context, "_", experiment, "_removing_outliers",".png")
      ggsave(filename, combined_pca_plot, width = 18, height = 12, dpi = 300, units = "cm")
      
    }, error = function(e) {
      message("Error plotting PCA for ", experiment, ", ", context, ": ", e$message)
    })
      
      
    } else if (experiment == "exp2"){
      tryCatch({
        # Perform PCA Analysis
        pca_results <- PCASamples(act_data_normed_fu, comp = c(1, 2), obj.return = T)
        pca_data <- as.data.frame(pca_results$x)
        pca_data$Sample <- rownames(pca_data)          
        pca_data$Sample <- as.factor(pca_data$Sample)
        
        # Calculate variance explained
        sdev <- pca_results$sdev
        variance <- sdev^2
        total_variance <- sum(variance)
        explained_variance <- variance / total_variance
        
        # Join with metadata
        pca_data <- pca_data %>%
          left_join(sample_metadata, by = c("Sample" = "SAMPLE_ID"))
        
        # Separate groups
        data_group_control <- pca_data %>% filter(Group == "B") # change this for exp2
        data_group_treatment <- pca_data %>% filter(Group == "C") # change this for exp2
        
        # Inner join to create *_control and *_treatment columns
        arrow_data <- inner_join(data_group_control, data_group_treatment, by = "Clone", suffix = c("_control", "_treatment")) %>%
          mutate(
            mid_x = (PC1_control + PC1_treatment) / 2,
            mid_y = (PC2_control + PC2_treatment) / 2
          )
        
        # Generate Hull data:
        df_hull <- pca_data %>%
          group_by(Group) %>%
          slice(chull(PC1, PC2))  # slice() keeps only hull points for each cluster
        
        
        # Generate PCA plot
        combined_pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
          geom_point(aes(color = Group), size = 5) +
          # geom_mark_hull(concavity = 5, expand = 0, radius = 0, aes(fill=Group, colour = Group)) + # old ggforce code
          geom_polygon(
            data = df_hull, 
            aes(x = PC1, y = PC2, fill = Group), 
            alpha = 0.2, 
            color = NA
          ) +
          geom_segment(
            data = arrow_data,
            aes(x = PC1_control, y = PC2_control, xend = PC1_treatment, yend = PC2_treatment),
            arrow = arrow(length = unit(0.2, "inches")),
            color = "black"
          ) +
          geom_label(data = arrow_data, aes(x = mid_x, y = mid_y, label = Clone), size = 3) +
          scale_color_manual(values = c("B" = "goldenrod2", "C" = "green")) +
          scale_fill_manual(values = c("B" = "goldenrod2", "C" = "green")) +
          theme_bw() +
          labs(
            title = "PCA Plot: PCs 1-2",
            x = paste0("PC", 1, " (", round(explained_variance[1] * 100, 1), "%)"),
            y = paste0("PC", 2, " (", round(explained_variance[2] * 100, 1), "%)"),
            color = "Group"
          )
        
        filename <- paste0("PCA_", 1, "-", 2, "_plot_", context, "_", experiment, "_removing_outliers",".png")
        ggsave(filename, combined_pca_plot, width = 18, height = 12, dpi = 300, units = "cm")
        
      }, error = function(e) {
        message("Error plotting PCA for ", experiment, ", ", context, ": ", e$message)
      })
      
    }
  
    ## Perform Differential Methylation Analyses ----
    # Differential methylation analyses
    dmb_data_exp <- NULL
    tryCatch({
      dmb_data_exp <- calculateDiffMeth(
        act_data_normed_fu,
        overdispersion = "MN",
        mc.cores = 32,
        suffix = paste0(context, "fu3", "odMNtestC")
      )
    }, error = function(e) {
      message("Error in differential methylation calculation for ", experiment, ", ", context, ", ", region, ": ", e$message)
    })
    
    if (is.null(dmb_data_exp)) {
      rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
      gc()
      next
    }
    
    ## Differential Methylation Data manipuatlion ----
    diffMethData <- NULL
    act_diff_data <- NULL
    act_diff_data_df_ann <- NULL
    tryCatch({
      
      # Get the data with the thresholds chosen
      diffMethData <- getMethylDiff(dmb_data_exp, 
                                    qvalue = 0.99, 
                                    difference = 5, 
                                    type = "all") # diff = 5 to prevent falling over.
      
      # change object to a df object to manipulate data
      act_diff_data <- getData(diffMethData)
      
      # !! APPLY THE P-ADJUSTMENT HERE !!
      act_diff_data$p_fdr = p.adjust(act_diff_data$pvalue, method = "fdr")
      
      # Load in matching file to be able to match to annotated genome file;
      methylation_matching <- read.table("../../Data/methylation_matching_file.tsv", 
                                         header = FALSE, 
                                         sep = "\t", 
                                         col.names = c("WHPX_id", "chr"))
      
      # left join the data files together from the act_diff_data:
      act_diff_data_matched <- left_join(act_diff_data, methylation_matching, by = "chr")

	  # save DMA data before joining to genome data:
	  write.table(
      act_diff_data_matched,
      file = paste0("Actinia_DMBs_d5_", experiment, "_", context, "_df_ann.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "")

      
      # change df to a gr object:
      act_diff_gr_matched <- as(act_diff_data_matched, "GRanges")
      
      # load in annotated genome
      genome <- read_gff3("/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Data2/combined_annotations.gff3")
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
      # 18. convert back to GRanges format
      genome <- as(genome, "GRanges")
      
      # left join the genome file to the data file:
      act_diff_gr_ann = join_overlap_left_directed(act_diff_gr_matched, genome)
      
      ## write table for data
      file.name = paste0("act_diff_df_ann", "_", experiment, "_", context, "_removing_outliers",".txt")
      act_diff_df_ann = as(act_diff_gr_ann, "data.frame")
      write.table(act_diff_df_ann, file.name, sep="\t", row.names=FALSE)
      
      # change back to gr object
      act_diff_gr_ann = as(act_diff_gr_ann, "GRanges")
      
    }, error = function(e) {
      message("Error in differential methylation steps for ", experiment, ", ", context, ": ", e$message)
    })
    
      # # Create volcano plot ----
    data <- act_diff_data
      data <- act_diff_df_ann
      data$diffmeth <- ifelse(data$meth.diff >= 0, "UP", "DOWN")
      data$diffmeth <- as.factor(data$diffmeth)
      
      volcano_plot <- ggplot(data = data, aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) +
        geom_point() +
        geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") +
        geom_vline(xintercept = c(5), col = "red", linetype = "dashed") +
        geom_hline(yintercept = c(1), col = "black", linetype = "dashed") +
        theme_classic() +
        theme(
          axis.title.y = element_text(
            face = "bold",
            margin = margin(0, 20, 0, 0),
            size = rel(1.1),
            color = "black"
          ),
          axis.title.x = element_text(
            hjust = 0.5,
            face = "bold",
            margin = margin(20, 0, 0, 0),
            size = rel(1.1),
            color = "black"
          ),
          plot.title = element_text(hjust = 0.5)
        ) +
        scale_color_manual(
          values = c("deepskyblue", "brown1"),
          labels = c("Down Methylated", "Up Methylated")
        ) +
        labs(
          x = "Differential methylation %",
          y = expression("-log"[10] * "p-adj")
        ) +
        ggtitle(paste("Differential Methylation of", context, "in", experiment, "removed outliers")) + 
        scale_x_continuous(limits = c(-40, 40), breaks = c(-40,-35,-30,-25,-20,-15, -10, -5, 0, 5, 10, 15, 20,25,30,35,40)) + 
        scale_y_continuous(limits = c(0, 3), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0,2.2,2.4,2.6,2.8,3.0))
      
      filename <- paste0("volcano_plot_", experiment, "_", context, "_removing_outliers",".png")
      ggsave(filename, volcano_plot, width = 18, height = 12, units = "cm")
      
    # Cleanup after each iteration
    rm(list = setdiff(ls(), c("experiments", "contexts",  "experiment", "context",  "base_path","genome","sample_metadata")))
    gc()
    
  }
  
  # Cleanup after each iteration
  rm(list = setdiff(ls(), c("experiments", "contexts","experiment", "context", "base_path","genome","sample_metadata")))
  gc()
  
} 
# Cleanup after each iteration
rm(list = setdiff(ls(), c("experiments", "contexts", "experiment", "context", "base_path","genome","sample_metadata")))
gc()

