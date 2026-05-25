#!/usr/bin/env Rscript
#SBATCH --job-name=regionalised_dma_no_outliers
#SBATCH --mem=100G
#SBATCH --partition=uoa-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

suppressPackageStartupMessages({
  library(methylKit)
  library(plyranges)
  library(GenomicRanges)
  library(IRanges)
  library(qqman)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# -------------------------------
# Project paths
# -------------------------------

project_dir <- "/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2"

regionalised_dir <- file.path(
  project_dir,
  "results",
  "regionalised_methylation",
  "aggregated_counts"
)

external_dir <- file.path(project_dir, "data", "external")
annotation_gff3 <- file.path(external_dir, "combined_annotations.gff3")

if (!file.exists(annotation_gff3)) {
  annotation_gff3_alt <- file.path(external_dir, "combined_annotation.gff3")
  if (file.exists(annotation_gff3_alt)) {
    annotation_gff3 <- annotation_gff3_alt
  }
}

metadata_file <- file.path(project_dir, "data", "metadata", "Methyl_sample_groups.txt")

results_dir <- file.path(project_dir, "results")

figure_dir <- file.path(results_dir, "figures", "regionalised_dma_no_outlier_pairs")
qc_figure_dir <- file.path(figure_dir, "methylkit_qc")
pca_figure_dir <- file.path(figure_dir, "pca")
volcano_figure_dir <- file.path(figure_dir, "volcano")

loadings_dir <- file.path(
  results_dir,
  "tables",
  "regionalised_dma_no_outlier_pairs",
  "pca_loadings"
)

model_output_dir <- file.path(
  results_dir,
  "model_outputs",
  "regionalised_dma_no_outlier_pairs"
)

tmp_methylkit_dir <- file.path(
  results_dir,
  "tmp",
  "regionalised_dma_no_outlier_pairs_methylkit_inputs"
)

dir.create(qc_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(loadings_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_methylkit_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "24"))

# -------------------------------
# Analysis settings
# -------------------------------

experiments <- c("exp1", "exp2")
contexts <- c("CpG", "CHG", "CHH")

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

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
# Load combined annotation
# -------------------------------

if (!file.exists(annotation_gff3)) {
  stop("Combined annotation GFF3 not found: ", annotation_gff3)
}

message("Loading combined annotation GFF3: ", annotation_gff3)

genome <- read_gff3(annotation_gff3)
genome_df <- as.data.frame(genome)

genome_df[] <- lapply(genome_df, function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  x[is.na(x)] <- "NA"
  x
})

if ("width" %in% names(genome_df)) {
  genome_df$width <- NULL
}

genome <- makeGRangesFromDataFrame(
  genome_df,
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

regions <- unique(as.character(mcols(genome)$type))
regions <- regions[!is.na(regions) & regions != "" & regions != "NA"]

preferred_regions <- c(
  "gene",
  "exon",
  "five_prime_UTR",
  "three_prime_UTR",
  "promoter",
  "downstream_region",
  "dispersed_repeat",
  "ncRNA_gene"
)

regions <- c(
  intersect(preferred_regions, regions),
  setdiff(regions, preferred_regions)
)

message("Detected regional feature types:")
message(paste(regions, collapse = ", "))

# -------------------------------
# Sample definitions: OUTLIER PAIRS EXCLUDED
# -------------------------------

get_sample_info <- function(experiment) {
  if (experiment == "exp1") {

    message("Excluding exp1 outlier pair: Sample_3-3_D and Sample_17-17_D")

    data.frame(
      sample_dir = c(
        "Sample_9-9_D",
        "Sample_13-13_D",
        "Sample_15-15_D",
        "Sample_16-16_D",
        "Sample_22-22_D",
        "Sample_27-27_D"
      ),
      file_sample_id = c("9", "13", "15", "16", "22", "27"),
      methylkit_sample_id = c("9", "13", "15", "16", "22", "27"),
      treatment = c(0, 0, 1, 0, 1, 1),
      ind_id = factor(c("15", "9", "15", "2", "9", "2")),
      stringsAsFactors = FALSE
    )

  } else if (experiment == "exp2") {

    message("Excluding exp2 outlier pair: Sample_8-8_D and Sample_32-18_redo_D")

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
      file_sample_id = c("6", "7", "10", "14", "20", "21", "29_11", "36_2"),
      methylkit_sample_id = c("6", "7", "10", "14", "20", "21", "11", "2"),
      treatment = c(0, 0, 0, 1, 1, 1, 0, 1),
      ind_id = factor(c("19", "8", "6", "8", "19", "12", "12", "6")),
      stringsAsFactors = FALSE
    )

  } else {
    stop("Unknown experiment: ", experiment)
  }
}

# -------------------------------
# Input preparation
# -------------------------------
# The regionalisation script writes files with metadata columns.
# methylKit expects bismarkCoverage-like six-column files:
# chr, start, end, percent_methylation, numCs, numTs.
# This function converts regionalised count files into temporary
# methylKit-compatible input files.

prepare_methylkit_region_file <- function(input_file, output_file) {
  dt <- fread(input_file)

  required_cols <- c("chr", "start", "end", "numCs", "numTs", "coverage")

  missing_cols <- setdiff(required_cols, names(dt))

  if (length(missing_cols) > 0) {
    stop(
      "Input file is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nFile: ",
      input_file
    )
  }

  dt <- dt[
    !is.na(chr) &
      !is.na(start) &
      !is.na(end) &
      !is.na(numCs) &
      !is.na(numTs) &
      !is.na(coverage)
  ]

  # methylKit cannot use regions with zero total coverage.
  dt <- dt[coverage > 0]

  if (nrow(dt) == 0) {
    stop("No covered regions remain after filtering: ", input_file)
  }

  dt[, percent_methylation := (numCs / coverage) * 100]

  out_dt <- dt[, .(
    chr,
    start,
    end,
    percent_methylation,
    numCs,
    numTs
  )]

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  fwrite(
    out_dt,
    file = output_file,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE,
    na = "NA"
  )

  output_file
}

build_file_list <- function(experiment, context, region, sample_info) {
  safe_region <- safe_name(region)

  input_files <- file.path(
    regionalised_dir,
    experiment,
    context,
    safe_region,
    paste0(
      sample_info$file_sample_id,
      "_",
      safe_region,
      "_aggregated_methylation_counts.txt"
    )
  )

  tmp_files <- file.path(
    tmp_methylkit_dir,
    experiment,
    context,
    safe_region,
    paste0(
      sample_info$file_sample_id,
      "_",
      safe_region,
      "_methylkit_input.txt"
    )
  )

  list(
    input_files = input_files,
    tmp_files = tmp_files
  )
}

# -------------------------------
# PCA plotting helper
# -------------------------------

plot_regional_pca <- function(act_data_normed_fu, experiment, context, region, sample_metadata) {
  pca_results <- PCASamples(
    act_data_normed_fu,
    comp = c(1, 2),
    obj.return = TRUE
  )

  loadings_abs <- abs(pca_results$rotation)
  safe_region <- safe_name(region)

  loadings_file <- file.path(
    loadings_dir,
    experiment,
    context,
    safe_region,
    paste0(
      "PCA_1-2_loadings_",
      context,
      "_",
      safe_region,
      "_",
      experiment,
      "_no_outlier_pairs.txt"
    )
  )

  dir.create(dirname(loadings_file), recursive = TRUE, showWarnings = FALSE)

  write.table(
    loadings_abs,
    loadings_file,
    sep = "\t",
    row.names = TRUE,
    quote = FALSE
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
    filter(n() >= 3) %>%
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
      title = paste("Regional PCA:", context, region, experiment, "outlier pairs excluded"),
      x = paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)"),
      color = "Group",
      fill = "Group"
    )

  pca_file <- file.path(
    pca_figure_dir,
    experiment,
    context,
    safe_region,
    paste0(
      "PCA_1-2_plot_",
      context,
      "_",
      safe_region,
      "_",
      experiment,
      "_no_outlier_pairs.png"
    )
  )

  dir.create(dirname(pca_file), recursive = TRUE, showWarnings = FALSE)

  ggsave(
    pca_file,
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
  sample_info <- get_sample_info(experiment)

  for (context in contexts) {
    for (region in regions) {
      safe_region <- safe_name(region)

      message("------------------------------------------------------------")
      message("Starting regionalised DMR analysis")
      message("Experiment: ", experiment)
      message("Context:    ", context)
      message("Region:     ", region)
      message("Outliers:   excluded")
      message("------------------------------------------------------------")

      paths <- build_file_list(
        experiment = experiment,
        context = context,
        region = region,
        sample_info = sample_info
      )

      missing_files <- paths$input_files[!file.exists(paths$input_files)]

      if (length(missing_files) > 0) {
        message("Missing regionalised input files:")
        message(paste(missing_files, collapse = "\n"))
        message("Skipping ", experiment, " / ", context, " / ", region)
        next
      }

      # Convert regionalised count files to methylKit-compatible six-column files
      converted_files <- character(length(paths$input_files))
      conversion_failed <- FALSE

      for (i in seq_along(paths$input_files)) {
        converted_files[i] <- tryCatch(
          prepare_methylkit_region_file(
            input_file = paths$input_files[i],
            output_file = paths$tmp_files[i]
          ),
          error = function(e) {
            message("Error preparing methylKit input:")
            message(e$message)
            conversion_failed <<- TRUE
            NA_character_
          }
        )
      }

      if (conversion_failed || any(is.na(converted_files))) {
        message(
          "Skipping ",
          experiment,
          " / ",
          context,
          " / ",
          region,
          " because input conversion failed."
        )
        next
      }

      file.list <- as.list(converted_files)
      sample.id <- as.list(sample_info$methylkit_sample_id)
      treatment <- sample_info$treatment
      covariates_df <- data.frame(ind_id = factor(sample_info$ind_id))

      # -------------------------------
      # Read methylation data
      # -------------------------------

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
        message(
          "Error reading methylation data for ",
          experiment,
          " / ",
          context,
          " / ",
          region
        )
        message("Error message: ", e$message)
      })

      if (is.null(act_data)) {
        gc()
        next
      }

      # -------------------------------
      # QC plots
      # -------------------------------

      tryCatch({
        qc_file <- file.path(
          qc_figure_dir,
          experiment,
          context,
          safe_region,
          paste0(
            "raw_methylation_plots_act_",
            context,
            "_",
            safe_region,
            "_",
            experiment,
            "_no_outlier_pairs.png"
          )
        )

        dir.create(dirname(qc_file), recursive = TRUE, showWarnings = FALSE)

        png(qc_file, width = 12, height = 8, units = "in", res = 300)
        lapply(act_data, getMethylationStats, plot = TRUE)
        dev.off()
      }, error = function(e) {
        message("Error plotting methylation stats: ", e$message)
        try(dev.off(), silent = TRUE)
      })

      tryCatch({
        coverage_file <- file.path(
          qc_figure_dir,
          experiment,
          context,
          safe_region,
          paste0(
            "raw_coverage_plots_act_",
            context,
            "_",
            safe_region,
            "_",
            experiment,
            "_no_outlier_pairs.png"
          )
        )

        dir.create(dirname(coverage_file), recursive = TRUE, showWarnings = FALSE)

        png(coverage_file, width = 20, height = 10, units = "in", res = 300)
        lapply(act_data, getCoverageStats, plot = TRUE)
        dev.off()
      }, error = function(e) {
        message("Error plotting coverage stats: ", e$message)
        try(dev.off(), silent = TRUE)
      })

      # -------------------------------
      # Filter and normalise
      # -------------------------------

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
        message("Error filtering coverage: ", e$message)
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
        message("Error normalising coverage: ", e$message)
      })

      if (is.null(act_data_f_norm)) {
        gc()
        next
      }

      # -------------------------------
      # Unite regional methylation data
      # -------------------------------

      options(datatable.allow.cartesian = TRUE)

      act_data_normed_fu <- NULL

      tryCatch({
        act_data_normed_fu <- methylKit::unite(
          act_data_f_norm,
          destrand = FALSE,
          min.per.group = 3L,
          suffix = paste0(context, "_normed_fu3")
        )
      }, error = function(e) {
        message("Error uniting methylation data: ", e$message)
      })

      if (is.null(act_data_normed_fu)) {
        gc()
        next
      }

      # -------------------------------
      # PCA
      # -------------------------------

      tryCatch({
        plot_regional_pca(
          act_data_normed_fu = act_data_normed_fu,
          experiment = experiment,
          context = context,
          region = region,
          sample_metadata = sample_metadata
        )
      }, error = function(e) {
        message("Error plotting PCA: ", e$message)
      })

      # -------------------------------
      # Differential regional methylation
      # -------------------------------

      dmb_data_exp <- NULL

      tryCatch({
        dmb_data_exp <- calculateDiffMeth(
          act_data_normed_fu,
          covariates = covariates_df,
          overdispersion = "MN",
          mc.cores = n_cores,
          suffix = paste0(context, "_fu3_", safe_region, "_odMNtestC_no_outlier_pairs")
        )
      }, error = function(e) {
        message("Error in differential methylation calculation: ", e$message)
      })

      if (is.null(dmb_data_exp)) {
        gc()
        next
      }

      # -------------------------------
      # Extract, annotate, plot, save
      # -------------------------------

      tryCatch({
        diffMethData <- getMethylDiff(
          dmb_data_exp,
          qvalue = 0.99,
          difference = 0,
          type = "all"
        )

        act_diff_data <- getData(diffMethData)

        act_diff_data$p_fdr <- p.adjust(
          act_diff_data$pvalue,
          method = "fdr"
        )

        act_diff_gr <- makeGRangesFromDataFrame(
          act_diff_data,
          keep.extra.columns = TRUE,
          ignore.strand = TRUE,
          seqnames.field = "chr",
          start.field = "start",
          end.field = "end"
        )

        genome_region <- genome[mcols(genome)$type == region]

        act_diff_gr_ann <- join_overlap_left_directed(
          act_diff_gr,
          genome_region
        )

        act_diff_data_df_ann <- as(act_diff_gr_ann, "data.frame")

        # -------------------------------
        # Volcano plot
        # -------------------------------

        volcano_data <- act_diff_data_df_ann

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
              region,
              "regions of",
              experiment,
              "with outlier pairs excluded"
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

        volcano_file <- file.path(
          volcano_figure_dir,
          experiment,
          context,
          safe_region,
          paste0(
            "volcano_plot_",
            experiment,
            "_",
            context,
            "_",
            safe_region,
            "_no_outlier_pairs_regions.png"
          )
        )

        dir.create(dirname(volcano_file), recursive = TRUE, showWarnings = FALSE)

        ggsave(
          volcano_file,
          volcano_plot,
          width = 18,
          height = 12,
          units = "cm",
          dpi = 300
        )

        # -------------------------------
        # Save outputs
        # -------------------------------

        output_dir <- file.path(
          model_output_dir,
          experiment,
          context,
          safe_region
        )

        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        rdata_file <- file.path(
          output_dir,
          paste0(
            "Actinia_DMRs_d5_",
            experiment,
            "_",
            context,
            "_",
            safe_region,
            "_no_outlier_pairs_df_ann.RData"
          )
        )

        txt_file <- file.path(
          output_dir,
          paste0(
            "Actinia_DMRs_d5_",
            experiment,
            "_",
            context,
            "_",
            safe_region,
            "_no_outlier_pairs_df_ann.txt"
          )
        )

        save(
          act_diff_data_df_ann,
          file = rdata_file
        )

        write.table(
          act_diff_data_df_ann,
          file = txt_file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE,
          na = ""
        )

        message("Saved regionalised DMR table: ", txt_file)

      }, error = function(e) {
        message(
          "Error in differential methylation processing for ",
          experiment,
          " / ",
          context,
          " / ",
          region,
          ": ",
          e$message
        )
      })

      rm(
        act_data,
        act_data_f,
        act_data_f_norm,
        act_data_normed_fu,
        dmb_data_exp
      )

      gc()

      message("Completed ", experiment, " / ", context, " / ", region)
    }
  }
}

message("All regionalised DMR analyses with outlier pairs excluded are complete.")