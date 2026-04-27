#!/usr/bin/env Rscript
#SBATCH --job-name=methylation_analysis
#SBATCH --mem=400G
#SBATCH --partition=uoa-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=21-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

suppressPackageStartupMessages({
  library(methylKit)
  library(plyranges)
  library(GenomicRanges)
  library(qqman)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# -------------------------------
# Project paths
# -------------------------------

project_dir <- "/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2"

bedgraph_dir <- file.path(project_dir, "data", "processed", "bedgraph")
metadata_file <- file.path(project_dir, "data", "metadata", "Methyl_sample_groups.txt")

reference_dir <- file.path(project_dir, "data", "reference")
annotation_gff3 <- file.path(reference_dir, "combined_annotations.gff3")
methylation_matching_file <- file.path(reference_dir, "methylation_matching_file.tsv")

results_dir <- file.path(project_dir, "results")
figure_dir <- file.path(results_dir, "figures")
qc_figure_dir <- file.path(figure_dir, "methylkit_qc")
pca_figure_dir <- file.path(figure_dir, "pca")
volcano_figure_dir <- file.path(figure_dir, "volcano")
model_output_dir <- file.path(results_dir, "model_outputs", "base_level_dma")

dir.create(qc_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)

setwd(model_output_dir)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))

# -------------------------------
# Analysis settings
# -------------------------------

experiments <- c("exp1", "exp2")
contexts <- c("CpG", "CHG", "CHH")

# -------------------------------
# Metadata
# -------------------------------

sample_metadata <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

sample_metadata <- sample_metadata %>%
  mutate(
    SAMPLE_ID = Label %>%
      sub("_D$", "", .) %>%
      sub("_redo$", "", .),
    SAMPLE_ID = as.factor(SAMPLE_ID),
    Group = as.factor(Group),
    Clone = as.factor(Clone)
  )

# -------------------------------
# Optional annotation inputs
# -------------------------------

genome <- NULL

if (file.exists(annotation_gff3)) {
  genome <- read_gff3(annotation_gff3)

  genome <- as.data.frame(genome)

  genome[] <- lapply(genome, function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }
    x[is.na(x)] <- "NA"
    x
  })

  genome <- as(genome, "GRanges")
} else {
  warning("Annotation GFF3 not found: ", annotation_gff3)
  warning("Differential methylation output will be written without GFF3 overlap annotation.")
}

methylation_matching <- NULL

if (file.exists(methylation_matching_file)) {
  methylation_matching <- read.table(
    methylation_matching_file,
    header = FALSE,
    sep = "\t",
    col.names = c("WHPX_id", "chr"),
    stringsAsFactors = FALSE
  )
} else {
  warning("Methylation matching file not found: ", methylation_matching_file)
  warning("Differential methylation output will be written without chromosome-name matching.")
}

# -------------------------------
# Sample definitions
# -------------------------------

get_sample_info <- function(experiment) {
  if (experiment == "exp1") {
    data.frame(
      sample_dir = c(
        "Sample_9-9_D",
        "Sample_13-13_D",
        "Sample_15-15_D",
        "Sample_16-16_D",
        "Sample_22-22_D",
        "Sample_27-27_D"
      ),
      sample_id = c("9", "13", "15", "16", "22", "27"),
      treatment = c(0, 0, 1, 0, 1, 1),
      ind_id = factor(c("15", "9", "15", "2", "9", "2")),
      stringsAsFactors = FALSE
    )
  } else if (experiment == "exp2") {
    data.frame(
      sample_dir = c(
        "Sample_6-6_D",
        "Sample_7-7_D",
        "Sample_10-10_D",
        "Sample_14-14_D",
        "Sample_20-20_D",
        "Sample_21-21_D",
        "Sample_29-11_redo_D",
        "Sample_36-2_redo_D"
      ),
      sample_id = c("6", "7", "10", "14", "20", "21", "11", "2"),
      treatment = c(0, 0, 0, 1, 1, 1, 0, 1),
      ind_id = factor(c("19", "8", "6", "8", "19", "12", "12", "6")),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown experiment: ", experiment)
  }
}

build_file_list <- function(sample_info, context) {
  file.path(
    bedgraph_dir,
    sample_info$sample_dir,
    paste0("meth_", context, "_cov_reads")
  )
}

# -------------------------------
# PCA plotting helper
# -------------------------------

plot_pca <- function(act_data_normed_fu, experiment, context, sample_metadata) {
  pca_results <- PCASamples(
    act_data_normed_fu,
    comp = c(1, 2),
    obj.return = TRUE
  )

  pca_data <- as.data.frame(pca_results$x)
  pca_data$Sample <- rownames(pca_data)
  pca_data$Sample <- as.factor(pca_data$Sample)

  sdev <- pca_results$sdev
  variance <- sdev^2
  explained_variance <- variance / sum(variance)

  pca_data <- pca_data %>%
    left_join(sample_metadata, by = c("Sample" = "SAMPLE_ID"))

  if (experiment == "exp1") {
    control_group <- "X"
    treatment_group <- "A"
    colour_values <- c("A" = "red", "X" = "blue")
  } else {
    control_group <- "B"
    treatment_group <- "C"
    colour_values <- c("B" = "goldenrod2", "C" = "green")
  }

  data_group_control <- pca_data %>% filter(Group == control_group)
  data_group_treatment <- pca_data %>% filter(Group == treatment_group)

  arrow_data <- inner_join(
    data_group_control,
    data_group_treatment,
    by = "Clone",
    suffix = c("_control", "_treatment")
  ) %>%
    mutate(
      mid_x = (PC1_control + PC1_treatment) / 2,
      mid_y = (PC2_control + PC2_treatment) / 2
    )

  df_hull <- pca_data %>%
    filter(!is.na(Group)) %>%
    group_by(Group) %>%
    slice(chull(PC1, PC2))

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
      aes(
        x = PC1_control,
        y = PC2_control,
        xend = PC1_treatment,
        yend = PC2_treatment
      ),
      arrow = arrow(length = unit(0.2, "inches")),
      color = "black"
    ) +
    geom_label(
      data = arrow_data,
      aes(x = mid_x, y = mid_y, label = Clone),
      size = 3
    ) +
    scale_color_manual(values = colour_values) +
    scale_fill_manual(values = colour_values) +
    theme_bw() +
    labs(
      title = paste("PCA plot:", context, experiment, "removing outliers"),
      x = paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)"),
      color = "Group",
      fill = "Group"
    )

  filename <- file.path(
    pca_figure_dir,
    paste0("PCA_1-2_plot_", context, "_", experiment, "_removing_outliers.png")
  )

  ggsave(
    filename,
    combined_pca_plot,
    width = 18,
    height = 12,
    dpi = 300,
    units = "cm"
  )
}

# -------------------------------
# Main loop
# -------------------------------

for (experiment in experiments) {
  for (context in contexts) {

    message("------------------------------------------------------------")
    message("Starting analysis: ", experiment, " / ", context)
    message("------------------------------------------------------------")

    sample_info <- get_sample_info(experiment)

    file.list <- build_file_list(sample_info, context)
    sample.id <- as.list(sample_info$sample_id)
    treatment <- sample_info$treatment
    covariates_df <- data.frame(ind_id = sample_info$ind_id)

    missing_files <- file.list[!file.exists(file.list)]

    if (length(missing_files) > 0) {
      message("Missing methylKit input files for ", experiment, " / ", context, ":")
      message(paste(missing_files, collapse = "\n"))
      message("Skipping this experiment/context.")
      next
    }

    # -------------------------------
    # Read methylation data
    # -------------------------------

    act_data <- NULL

    tryCatch({
      act_data <- methRead(
        as.list(file.list),
        sample.id = sample.id,
        assembly = "actinia",
        pipeline = "bismarkCoverage",
        treatment = treatment,
        context = context,
        dbtype = "tabix",
        mincov = 0
      )
    }, error = function(e) {
      message("Error reading methylation data for ", experiment, " / ", context)
      message("Error message: ", e$message)
    })

    if (is.null(act_data)) {
      gc()
      next
    }

    # -------------------------------
    # Methylation and coverage QC plots
    # -------------------------------

    tryCatch({
      filename <- file.path(
        qc_figure_dir,
        paste0("raw_methylation_plots_act_", context, "_", experiment, "_removing_outliers.png")
      )

      png(filename, width = 12, height = 8, units = "in", res = 300)
      lapply(act_data, getMethylationStats, plot = TRUE)
      dev.off()
    }, error = function(e) {
      message("Error plotting methylation stats for ", experiment, " / ", context, ": ", e$message)
      try(dev.off(), silent = TRUE)
    })

    tryCatch({
      filename <- file.path(
        qc_figure_dir,
        paste0("raw_coverage_plots_act_", context, "_", experiment, "_removing_outliers.png")
      )

      png(filename, width = 20, height = 10, units = "in", res = 300)
      lapply(act_data, getCoverageStats, plot = TRUE)
      dev.off()
    }, error = function(e) {
      message("Error plotting coverage stats for ", experiment, " / ", context, ": ", e$message)
      try(dev.off(), silent = TRUE)
    })

    # -------------------------------
    # Coverage filtering and normalisation
    # -------------------------------

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
      message("Error filtering coverage for ", experiment, " / ", context, ": ", e$message)
    })

    if (is.null(act_data_f)) {
      gc()
      next
    }

    act_data_f_norm <- NULL

    tryCatch({
      act_data_f_norm <- normalizeCoverage(
        act_data_f,
        method = "median"
      )
    }, error = function(e) {
      message("Error normalising coverage for ", experiment, " / ", context, ": ", e$message)
    })

    if (is.null(act_data_f_norm)) {
      gc()
      next
    }

    # -------------------------------
    # Unite methylation data
    # -------------------------------

    options(datatable.allow.cartesian = FALSE)

    act_data_normed_fu <- NULL

    tryCatch({
      act_data_normed_fu <- methylKit::unite(
        act_data_f_norm,
        destrand = FALSE,
        min.per.group = 3L,
        suffix = paste0(context, "_normed_fu3")
      )
    }, error = function(e) {
      message("Error uniting methylation data for ", experiment, " / ", context, ": ", e$message)
    })

    if (is.null(act_data_normed_fu)) {
      gc()
      next
    }

    # -------------------------------
    # Clustering plot
    # -------------------------------

    tryCatch({
      filename <- file.path(
        qc_figure_dir,
        paste0("clustering_plot_", context, "_", experiment, "_removing_outliers.png")
      )

      png(filename, width = 18, height = 12, units = "cm", res = 300)

      clusterSamples(
        act_data_normed_fu,
        filterByQuantile = FALSE,
        sd.threshold = 0.1,
        dist = "correlation",
        method = "ward.D",
        plot = TRUE
      )

      dev.off()
    }, error = function(e) {
      message("Error generating clustering plot for ", experiment, " / ", context, ": ", e$message)
      try(dev.off(), silent = TRUE)
    })

    # -------------------------------
    # PCA plot
    # -------------------------------

    tryCatch({
      plot_pca(
        act_data_normed_fu = act_data_normed_fu,
        experiment = experiment,
        context = context,
        sample_metadata = sample_metadata
      )
    }, error = function(e) {
      message("Error plotting PCA for ", experiment, " / ", context, ": ", e$message)
    })

    # -------------------------------
    # Differential methylation analysis
    # -------------------------------

    dmb_data_exp <- NULL

    tryCatch({
      dmb_data_exp <- calculateDiffMeth(
        act_data_normed_fu,
        covariates = covariates_df,
        overdispersion = "MN",
        mc.cores = n_cores,
        suffix = paste0(context, "_fu3_odMNtestC")
      )
    }, error = function(e) {
      message("Error in differential methylation calculation for ", experiment, " / ", context, ": ", e$message)
    })

    if (is.null(dmb_data_exp)) {
      gc()
      next
    }

    # -------------------------------
    # Extract, adjust, match, and annotate DMB data
    # -------------------------------

    act_diff_data <- NULL
    act_diff_df_ann <- NULL

    tryCatch({
      diffMethData <- getMethylDiff(
        dmb_data_exp,
        qvalue = 0.99,
        difference = 5,
        type = "all"
      )

      act_diff_data <- getData(diffMethData)

      act_diff_data$p_fdr <- p.adjust(
        act_diff_data$pvalue,
        method = "fdr"
      )

      if (!is.null(methylation_matching)) {
        act_diff_data_matched <- left_join(
          act_diff_data,
          methylation_matching,
          by = "chr"
        )
      } else {
        act_diff_data_matched <- act_diff_data
      }

      raw_output_file <- file.path(
        model_output_dir,
        paste0("Actinia_DMBs_d5_", experiment, "_", context, "_removing_outliers_raw.txt")
      )

      write.table(
        act_diff_data_matched,
        file = raw_output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
      )

      if (!is.null(genome)) {
        act_diff_gr_matched <- as(act_diff_data_matched, "GRanges")

        act_diff_gr_ann <- join_overlap_left_directed(
          act_diff_gr_matched,
          genome
        )

        act_diff_df_ann <- as(act_diff_gr_ann, "data.frame")
      } else {
        act_diff_df_ann <- act_diff_data_matched
      }

      annotated_output_file <- file.path(
        model_output_dir,
        paste0("act_diff_df_ann_", experiment, "_", context, "_removing_outliers.txt")
      )

      write.table(
        act_diff_df_ann,
        annotated_output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
      )

    }, error = function(e) {
      message("Error in differential methylation processing for ", experiment, " / ", context, ": ", e$message)
    })

    if (is.null(act_diff_data) || is.null(act_diff_df_ann)) {
      gc()
      next
    }

    # -------------------------------
    # Volcano plot
    # -------------------------------

    tryCatch({
      volcano_data <- act_diff_df_ann

      volcano_data$diffmeth <- ifelse(
        volcano_data$meth.diff >= 0,
        "UP",
        "DOWN"
      )

      volcano_data$diffmeth <- as.factor(volcano_data$diffmeth)

      volcano_plot <- ggplot(
        data = volcano_data,
        aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)
      ) +
        geom_point() +
        geom_vline(xintercept = -5, col = "blue", linetype = "dashed") +
        geom_vline(xintercept = 5, col = "red", linetype = "dashed") +
        geom_hline(yintercept = 1, col = "black", linetype = "dashed") +
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
          labels = c("Down methylated", "Up methylated")
        ) +
        labs(
          x = "Differential methylation %",
          y = expression("-log"[10] * " FDR-adjusted p-value"),
          col = "Direction"
        ) +
        ggtitle(
          paste(
            "Differential methylation of",
            context,
            "in",
            experiment,
            "with outliers removed"
          )
        ) +
        scale_x_continuous(
          limits = c(-40, 40),
          breaks = seq(-40, 40, by = 5)
        ) +
        scale_y_continuous(
          limits = c(0, 3),
          breaks = seq(0, 3, by = 0.2)
        )

      filename <- file.path(
        volcano_figure_dir,
        paste0("volcano_plot_", experiment, "_", context, "_removing_outliers.png")
      )

      ggsave(
        filename,
        volcano_plot,
        width = 18,
        height = 12,
        units = "cm",
        dpi = 300
      )

    }, error = function(e) {
      message("Error creating volcano plot for ", experiment, " / ", context, ": ", e$message)
    })

    rm(
      act_data,
      act_data_f,
      act_data_f_norm,
      act_data_normed_fu,
      dmb_data_exp,
      act_diff_data,
      act_diff_df_ann
    )

    gc()

    message("Completed analysis: ", experiment, " / ", context)
  }
}

message("All methylation analyses complete.")