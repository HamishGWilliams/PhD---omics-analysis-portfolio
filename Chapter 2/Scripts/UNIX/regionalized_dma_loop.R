#!/usr/bin/env Rscript
#SBATCH --job-name=DMR_analysis
#SBATCH --mem=400G
#SBATCH --partition=uoa-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
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
library(qqman) # for Q-Q and Manhattan plots  
library(dplyr) # for data manipulation  
library(ggplot2) # for plotting  
library(data.table) # for data table manipulation

# Set base_path
base_path <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results"
# change to base_path
setwd(base_path)

# Define Parameters for looping:
experiments <- c("exp1","exp2")
contexts <- c("CpG","CHG","CHH") 

	# Running CHH context with 200G to try and brute force the analyses (28th April 2025)
	# Doing this because I think there might be some interesting observations to be had with these cytosines

	# contexts <- c("CHH")

# load genome
genome <- read_gff3("../../Data2/combined_annotations.gff3")

# extract region unique names for looping
regions <- unique(genome$type)

# load sample_metadata for adding data to plots and analysis:
sample_metadata <- read.table("./sample_metadata.txt", header = T)
sample_metadata$SAMPLE_ID <- as.factor(sample_metadata$SAMPLE_ID)

# ----
# Start Loop

for (experiment in experiments) {
  for (context in contexts) {
    for (region in regions) {
      
      # Reactively make new path to specific folder:
      exp_context_region_path <- file.path(base_path, experiment, context, region)
      # set new directory
      setwd(exp_context_region_path)
      # Ensure objects are empty:
      file.list <- NULL
      sample.id <- NULL
      treatment <- NULL
      
      # Load sample files depending on experiment
      if (experiment == "exp1") {
        file.list = list(
          # paste0("./3_", region, "_aggregated_methylation_counts.txt"), # excluding outlier
          paste0("./9_", region, "_aggregated_methylation_counts.txt"),
          paste0("./13_", region, "_aggregated_methylation_counts.txt"),
          paste0("./15_", region, "_aggregated_methylation_counts.txt"),
          paste0("./16_", region, "_aggregated_methylation_counts.txt"),
          # paste0("./17_", region, "_aggregated_methylation_counts.txt"), # excluding outlier
          paste0("./22_", region, "_aggregated_methylation_counts.txt"),
          paste0("./27_", region, "_aggregated_methylation_counts.txt")
        )
        sample.id = list("9","13","15","16","22","27") # removed outliers from here too
        treatment = c(0, 0, 1, 0, 1, 1) # also removed outliers from here, otherwise the script wont run
      # Build covariate df for "fixed paired analysis" (changes made based on comments 14/12/25)
      covariates_exp1 = data.frame(ind_id = factor(c("15", "9", "15", "2", "9", "2")))

      } else if (experiment == "exp2") { # kind of redundant for exp1 script since exp2 isnt included in loop
        file.list = list(
          paste0("./6_", region, "_aggregated_methylation_counts.txt"),
          paste0("./7_", region, "_aggregated_methylation_counts.txt"),
          # paste0("./8_", region, "_aggregated_methylation_counts.txt"),
          paste0("./10_", region, "_aggregated_methylation_counts.txt"),
          paste0("./14_", region, "_aggregated_methylation_counts.txt"),
          paste0("./20_", region, "_aggregated_methylation_counts.txt"),
          paste0("./21_", region, "_aggregated_methylation_counts.txt"),
          paste0("./29_11_", region, "_aggregated_methylation_counts.txt"),
          # paste0("./32_18_", region, "_aggregated_methylation_counts.txt"),
          paste0("./36_2_", region, "_aggregated_methylation_counts.txt")
        )
        sample.id = list("6", "7","10", "14", "20", "21", "11", "2")
        treatment = c(0, 0, 0, 1, 1, 1, 0, 1)
      # Build covariate df for "fixed paired analysis" (changes made based on comments 14/12/25)
      covariates_exp2 = data.frame(ind_id = factor(c("19", "8", "6", "8", "19", "12", "12", "6")))

      }
      
      # Read the methylation data
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
          mincov = 0,
          resolution = "region"
        )
      }, error = function(e) {
        message("Error reading methylation data for experiment: ", experiment, ", context: ", context, ", region: ", region)
        message("Error message: ", e$message)
      })
      
      # If methRead failed, no point in continuing to next steps
      if (is.null(act_data)) {
        # Clean up and continue to next iteration
        rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
        gc()
        next
      }
      
      # Plot methylation statistics
      tryCatch({
        filename <- paste0("raw_methylation_plots_act_", context, "_", region, "_regions_", experiment, "without_outlier",".png")
        png(filename, width = 12, height = 8, units = "in", res = 300)
        lapply(act_data, getMethylationStats, plot = TRUE)
        dev.off()
      }, error = function(e) {
        message("Error plotting methylation stats for ", experiment, ", ", context, ", ", region, ": ", e$message)
        dev.off()
      })
      
      # Plot raw coverage statistics
      tryCatch({
        filename <- paste0("raw_coverage_plots_act_", context, "_", region, "_regions_", experiment, "without_outlier", ".png")
        png(filename, width = 20, height = 10, units = "in", res = 300)
        lapply(act_data, getCoverageStats, plot = TRUE)
        dev.off()
      }, error = function(e) {
        message("Error plotting coverage stats for ", experiment, ", ", context, ", ", region, ": ", e$message)
        dev.off()
      })
      
      # Filter for coverage
      act_data_f <- NULL
      tryCatch({
        act_data_f <- filterByCoverage(
          act_data,
          hi.count = 10000,
          hi.perc = NULL,
          lo.count = 5,
          lo.perc = NULL,
          suffix = paste0(context, "_f")
        )
      }, error = function(e) {
        message("Error filtering coverage for ", experiment, ", ", context, ", ", region, ": ", e$message)
      })
      
      if (is.null(act_data_f)) {
        rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
        gc()
        next
      }
      
      # Normalize coverage
      act_data_f_norm <- NULL
      tryCatch({
        act_data_f_norm <- normalizeCoverage(act_data_f, method = "median")
      }, error = function(e) {
        message("Error normalizing coverage for ", experiment, ", ", context, ", ", region, ": ", e$message)
      })
      
      if (is.null(act_data_f_norm)) {
        rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
        gc()
        next
      }
      
      # Set data.table option
      options(datatable.allow.cartesian = T) # necessary to generate results in the first place, however it does introduce some duplicate results in some instances.
      
      # Unite methylation results
      act_data_normed_fu <- NULL
      tryCatch({
        act_data_normed_fu <- methylKit::unite(
          act_data_f_norm,
          destrand = FALSE,
          min.per.group = 3L,
          suffix = paste0(context, "_normed_fu3")
        )
      }, error = function(e) {
        message("Error uniting methylation data for ", experiment, ", ", context, ", ", region, ": ", e$message)
      })
      
      if (is.null(act_data_normed_fu)) {
        rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
        gc()
        next
      }	
      
        # Perform PCA Analysis
        pca_results <- PCASamples(act_data_normed_fu, comp = c(1, 2), obj.return = T)
        
        # Catch the loading values for later inspection:
        loadings <- pca_results$rotation
        loadings_abs <- abs(loadings)
        
        filename <- paste0("PCA_1-2_", "loading_values_", context, "_", region, "_", experiment, ".txt")
        write.table(loadings_abs, filename, sep = "\t", row.names = T)
        
        # Capture loadings data and save
        loadings <- pca_results$rotation
        loadings_abs <- abs(loadings)
        filename = paste0("PCA_", 1, "-", 2, "_loadings_", context, "_", experiment, "_removed_outliers", ".txt")
        write.table(loadings_abs, filename, sep = "\t", row.names = TRUE)
        
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
        
        # Different loop depending on experiment:
        if (experiment == "exp1"){
          
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
            scale_color_manual(values = c("A" = "red", "X" = "blue")) +
            scale_fill_manual(values = c("A" = "red", "X" = "blue")) +
            theme_bw() +
            labs(
              title = "PCA Plot: PCs 1-2",
              x = paste0("PC", 1, " (", round(explained_variance[1] * 100, 1), "%)"),
              y = paste0("PC", 2, " (", round(explained_variance[2] * 100, 1), "%)"),
              color = "Group"
            )
          
          filename <- paste0("PCA_", 1, "-", 2, "_plot_", context, "_", region, "_", experiment, ".png")
          ggsave(filename, combined_pca_plot, width = 18, height = 12, dpi = 300, units = "cm")
          
        } else if (experiment == "exp2"){
          
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
          
          
          combined_pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
            geom_point(aes(color = Group), size = 5) +
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
          
          filename <- paste0("PCA_", 1, "-", 2, "_plot_", context, "_", region, "_", experiment, ".png")
          ggsave(filename, combined_pca_plot, width = 18, height = 12, dpi = 300, units = "cm")
        }
      
      # Differential methylation analyses
	cov_name <- paste0("covariates_", experiment)
      dmb_data_exp <- NULL
      tryCatch({
        dmb_data_exp <- calculateDiffMeth(
          act_data_normed_fu,
	  covariates = get(cov_name, envir = .GlobalEnv),
          overdispersion = "MN",
          mc.cores = 32,
          suffix = paste0(context, "fu3", region, "odMNtestC")
        )
      }, error = function(e) {
        message("Error in differential methylation calculation for ", experiment, ", ", context, ", ", region, ": ", e$message)
      })
      
      if (is.null(dmb_data_exp)) {
        rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
        gc()
        next
      }
      
	  ## Differntial Methylation Data manipulation
      diffMethData <- NULL
      act_diff_data <- NULL
      act_diff_data_df_ann <- NULL
      tryCatch({
        diffMethData <- getMethylDiff(dmb_data_exp, qvalue = 0.99, difference = 0, type = "all") # changed difference = 0 since we have power to compute all data per iteration
        
        act_diff_data <- getData(diffMethData)
          
        # !! IMPORTANT !! P VALUE ADJUSTMENTS ARE MADE HERE TO AVOIDE OVER INFLATION OF ROWS
        act_diff_data$p_fdr <- p.adjust(act_diff_data$pvalue, method = "fdr")
        
        act_diff_data_gr <- as(act_diff_data, "GRanges")
        
        genome_ann <- read_gff3("../../../../../Data2/combined_annotations.gff3")
        genome_ann <- genome_ann[genome_ann$type == region]
        act_diff_data_gr_ann <- join_overlap_left_directed(act_diff_data_gr, genome_ann)
        act_diff_data_df_ann <- as(act_diff_data_gr_ann, "data.frame")
        
        # Create volcano plot
        data <- act_diff_data_df_ann
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
          ggtitle(paste("Differential Methylation of", context, "in", region, "regions of", experiment, "without_outlier")) + 
          scale_x_continuous(limits = c(-40, 40), breaks = c(-40,-35,-30,-25,-20,-15, -10, -5, 0, 5, 10, 15, 20,25,30,35,40)) + 
          scale_y_continuous(limits = c(0, 3), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0,2.2,2.4,2.6,2.8,3.0))
        
        filename <- paste0("volcano_plot_", experiment, "_", context, "_", region, "without_outlier", "_regions.png")
        ggsave(filename, volcano_plot, width = 18, height = 12, units = "cm")
        
        # Save final annotated data frame
        save(act_diff_data_df_ann, file = file.path(exp_context_region_path, paste0("Actinia_DMBs_d5_", experiment, "_", context, "_", region, "without_outlier","_df_ann.RData")))
        
        write.table(
          act_diff_data_df_ann,
          file = paste0("Actinia_DMBs_d5_", experiment, "_", context, "_", region, "without_outlier","_df_ann.txt"),
          sep = "\t",
          quote = FALSE,
          row.names = FALSE,
          na = ""
        )
      }, error = function(e) {
        message("Error in differential methylation steps for ", experiment, ", ", context, ", ", region, ": ", e$message)
      })
      
      # Cleanup after each iteration
      rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
      gc()
      
    } # end region loop
    rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
    gc()
    
  } # end context loop
  rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
  gc()
  
} # end experiment loop
rm(list = setdiff(ls(), c("experiments", "contexts", "regions", "experiment", "context", "region", "base_path","genome","sample_metadata")))
gc()
