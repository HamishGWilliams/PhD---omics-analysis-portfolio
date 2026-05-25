#!/usr/bin/env Rscript
#SBATCH --job-name=base_level_dma_covariates
#SBATCH --mem=200G
#SBATCH --partition=uoa-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
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

external_dir <- file.path(project_dir, "data", "external")
annotation_gff3 <- file.path(external_dir, "combined_annotations.gff3")
methylation_matching_file <- file.path(external_dir, "methylation_matching_file.tsv")

# Fallback if the matching file is still stored in data/reference
if (!file.exists(methylation_matching_file)) {
  methylation_matching_file <- file.path(
    project_dir,
    "data",
    "reference",
    "methylation_matching_file.tsv"
  )
}

results_dir <- file.path(project_dir, "results")
figure_dir <- file.path(results_dir, "figures")

qc_figure_dir <- file.path(
  figure_dir,
  "methylkit_qc_base_level_covariates_no_outlier_pairs"
)

pca_figure_dir <- file.path(
  figure_dir,
  "pca_base_level_covariates_no_outlier_pairs"
)

volcano_figure_dir <- file.path(
  figure_dir,
  "volcano_base_level_covariates_no_outlier_pairs"
)

model_output_dir <- file.path(
  results_dir,
  "model_outputs",
  "base_level_dma_covariates_no_outlier_pairs"
)

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

# Normal context run
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
# Helper functions
# -------------------------------

clean_df_for_output <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  df[] <- lapply(df, function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }

    if (is.list(x)) {
      x <- vapply(
        x,
        function(y) {
          if (length(y) == 0 || all(is.na(y))) {
            NA_character_
          } else {
            paste(as.character(y), collapse = ";")
          }
        },
        character(1)
      )
    }

    x
  })

  df
}

add_annotation_chr_column <- function(diff_df, methylation_matching) {
  diff_df <- as.data.frame(diff_df, stringsAsFactors = FALSE)

  if (!"chr" %in% names(diff_df)) {
    stop("Differential methylation table does not contain a chr column.")
  }

  diff_df$chr <- as.character(diff_df$chr)
  diff_df$chr_for_annotation <- diff_df$chr

  if (is.null(methylation_matching)) {
    return(diff_df)
  }

  if (!all(c("WHPX_id", "chr") %in% names(methylation_matching))) {
    warning(
      "Methylation matching file must contain WHPX_id and chr columns. ",
      "Using original chr values."
    )
    return(diff_df)
  }

  mm <- methylation_matching %>%
    mutate(
      WHPX_id = as.character(WHPX_id),
      chr = as.character(chr)
    ) %>%
    distinct(WHPX_id, chr, .keep_all = TRUE)

  n_match_whpx <- sum(diff_df$chr %in% mm$WHPX_id, na.rm = TRUE)
  n_match_chr <- sum(diff_df$chr %in% mm$chr, na.rm = TRUE)

  if (n_match_whpx > 0 && n_match_whpx >= n_match_chr) {
    message(
      "Using methylation_matching_file.tsv to map methylKit chr values ",
      "from WHPX_id to annotation chr."
    )

    map <- mm %>%
      distinct(WHPX_id, .keep_all = TRUE) %>%
      transmute(
        chr = WHPX_id,
        chr_for_annotation_match = chr
      )

    diff_df <- left_join(diff_df, map, by = "chr")

    diff_df$chr_for_annotation <- ifelse(
      !is.na(diff_df$chr_for_annotation_match),
      diff_df$chr_for_annotation_match,
      diff_df$chr
    )

    diff_df$chr_for_annotation_match <- NULL

  } else if (n_match_chr > 0) {
    message(
      "methylKit chr values already match the chr column in ",
      "methylation_matching_file.tsv."
    )

    map <- mm %>%
      distinct(chr, .keep_all = TRUE) %>%
      transmute(
        chr = chr,
        WHPX_id_match = WHPX_id
      )

    diff_df <- left_join(diff_df, map, by = "chr")

  } else {
    warning(
      "No chromosome matches found in methylation_matching_file.tsv. ",
      "Using original chr values."
    )
  }

  diff_df
}

annotate_diff_table <- function(diff_df, genome, seq_col = "chr_for_annotation") {
  diff_df <- as.data.frame(diff_df, stringsAsFactors = FALSE)

  if (is.null(genome) || nrow(diff_df) == 0) {
    return(diff_df)
  }

  required_cols <- c(seq_col, "start", "end")
  missing_cols <- setdiff(required_cols, names(diff_df))

  if (length(missing_cols) > 0) {
    warning(
      "Cannot annotate differential methylation table. Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
    return(diff_df)
  }

  ann_df <- as.data.frame(genome)
  ann_df <- clean_df_for_output(ann_df)

  names(ann_df) <- paste0("annotation_", names(ann_df))
  names(ann_df) <- make.unique(names(ann_df), sep = "_")

  seqnames_vec <- as.character(diff_df[[seq_col]])
  start_vec <- suppressWarnings(as.integer(diff_df$start))
  end_vec <- suppressWarnings(as.integer(diff_df$end))

  valid <- !is.na(seqnames_vec) &
    !is.na(start_vec) &
    !is.na(end_vec) &
    start_vec <= end_vec

  if (!any(valid)) {
    warning("No valid genomic coordinates available for annotation.")

    empty_ann <- ann_df[rep(NA_integer_, nrow(diff_df)), , drop = FALSE]
    return(cbind(diff_df, empty_ann))
  }

  query_index <- which(valid)

  query_gr <- GenomicRanges::GRanges(
    seqnames = seqnames_vec[valid],
    ranges = IRanges::IRanges(
      start = start_vec[valid],
      end = end_vec[valid]
    ),
    strand = "*"
  )

  hits <- suppressWarnings(
    GenomicRanges::findOverlaps(
      query_gr,
      genome,
      ignore.strand = TRUE
    )
  )

  message("GFF3 overlap hits: ", length(hits))

  if (length(hits) == 0) {
    empty_ann <- ann_df[rep(NA_integer_, nrow(diff_df)), , drop = FALSE]
    return(cbind(diff_df, empty_ann))
  }

  hit_query_rows <- query_index[S4Vectors::queryHits(hits)]
  hit_subject_rows <- S4Vectors::subjectHits(hits)

  hit_out <- cbind(
    data.frame(.diff_original_row = hit_query_rows),
    diff_df[hit_query_rows, , drop = FALSE],
    ann_df[hit_subject_rows, , drop = FALSE]
  )

  unmatched_rows <- setdiff(seq_len(nrow(diff_df)), unique(hit_query_rows))

  if (length(unmatched_rows) > 0) {
    unmatched_out <- cbind(
      data.frame(.diff_original_row = unmatched_rows),
      diff_df[unmatched_rows, , drop = FALSE],
      ann_df[rep(NA_integer_, length(unmatched_rows)), , drop = FALSE]
    )

    out <- rbind(hit_out, unmatched_out)
  } else {
    out <- hit_out
  }

  out <- out[order(out$.diff_original_row), , drop = FALSE]
  out$.diff_original_row <- NULL
  rownames(out) <- NULL

  out
}

# -------------------------------
# Load combined annotation
# -------------------------------

genome <- NULL

if (file.exists(annotation_gff3)) {
  message("Loading combined annotation GFF3: ", annotation_gff3)

  genome <- read_gff3(annotation_gff3)

  genome_mcols <- clean_df_for_output(S4Vectors::mcols(genome))
  S4Vectors::mcols(genome) <- S4Vectors::DataFrame(genome_mcols)

} else {
  warning("Annotation GFF3 not found: ", annotation_gff3)
  warning(
    "Differential methylation output will be written without GFF3 ",
    "overlap annotation."
  )
}

# -------------------------------
# Load methylation chromosome matching file
# -------------------------------

methylation_matching <- NULL

if (file.exists(methylation_matching_file)) {
  message("Loading methylation matching file: ", methylation_matching_file)

  methylation_matching <- read.table(
    methylation_matching_file,
    header = FALSE,
    sep = "\t",
    col.names = c("WHPX_id", "chr"),
    stringsAsFactors = FALSE
  )
} else {
  warning("Methylation matching file not found: ", methylation_matching_file)
  warning(
    "Differential methylation output will be written without chromosome-name ",
    "matching."
  )
}

# -------------------------------
# Sample definitions
# Outlier pairs excluded
# -------------------------------

get_sample_info <- function(experiment) {
  if (experiment == "exp1") {
    message(
      "Using exp1 with outlier pair excluded: ",
      "Sample_3-3_D and Sample_17-17_D"
    )

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
    message(
      "Using exp2 with outlier pair excluded: ",
      "Sample_8-8_D and Sample_32-18_redo_D"
    )

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

  data_group_control <- pca_data %>%
    filter(Group == control_group)

  data_group_treatment <- pca_data %>%
    filter(Group == treatment_group)

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
      title = paste(
        "PCA plot:",
        context,
        experiment,
        "covariate model, outlier pairs excluded"
      ),
      x = paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)"),
      color = "Group",
      fill = "Group"
    )

  filename <- file.path(
    pca_figure_dir,
    paste0(
      "PCA_1-2_plot_",
      context,
      "_",
      experiment,
      "_covariates_no_outlier_pairs.png"
    )
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
    message("Starting base-level DMA with covariates")
    message("Experiment: ", experiment)
    message("Context:    ", context)
    message("Outliers:   excluded")
    message("------------------------------------------------------------")

    sample_info <- get_sample_info(experiment)

    file.list <- build_file_list(sample_info, context)
    sample.id <- as.list(sample_info$sample_id)
    treatment <- sample_info$treatment
    covariates_df <- data.frame(ind_id = factor(sample_info$ind_id))

    if (nrow(covariates_df) != length(file.list)) {
      stop(
        "Covariate dataframe rows do not match number of methylation files for ",
        experiment,
        " / ",
        context
      )
    }

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
    # QC plots
    # -------------------------------

    tryCatch({
      filename <- file.path(
        qc_figure_dir,
        paste0(
          "raw_methylation_plots_act_",
          context,
          "_",
          experiment,
          "_covariates_no_outlier_pairs.png"
        )
      )

      png(filename, width = 12, height = 8, units = "in", res = 300)
      lapply(act_data, getMethylationStats, plot = TRUE)
      dev.off()
    }, error = function(e) {
      message(
        "Error plotting methylation stats for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
      try(dev.off(), silent = TRUE)
    })

    tryCatch({
      filename <- file.path(
        qc_figure_dir,
        paste0(
          "raw_coverage_plots_act_",
          context,
          "_",
          experiment,
          "_covariates_no_outlier_pairs.png"
        )
      )

      png(filename, width = 20, height = 10, units = "in", res = 300)
      lapply(act_data, getCoverageStats, plot = TRUE)
      dev.off()
    }, error = function(e) {
      message(
        "Error plotting coverage stats for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
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
      message(
        "Error filtering coverage for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
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
      message(
        "Error normalising coverage for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
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
      message(
        "Error uniting methylation data for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
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
        paste0(
          "clustering_plot_",
          context,
          "_",
          experiment,
          "_covariates_no_outlier_pairs.png"
        )
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
      message(
        "Error generating clustering plot for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
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
      message(
        "Error plotting PCA for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
    })

    # -------------------------------
    # Differential methylation analysis
    # Covariate model
    # -------------------------------

    dmb_data_exp <- NULL

    tryCatch({
      dmb_data_exp <- calculateDiffMeth(
        act_data_normed_fu,
        covariates = covariates_df,
        overdispersion = "MN",
        mc.cores = n_cores,
        suffix = paste0(context, "_fu3_odMN_covariates_no_outlier_pairs")
      )
    }, error = function(e) {
      message(
        "Error in differential methylation calculation for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
    })

    if (is.null(dmb_data_exp)) {
      gc()
      next
    }

    # -------------------------------
    # Extract, adjust, match, annotate
    # -------------------------------

    act_diff_data <- NULL
    act_diff_data_matched <- NULL
    act_diff_df_ann <- NULL

    tryCatch({
      diffMethData <- getMethylDiff(
        dmb_data_exp,
        qvalue = 0.99,
        difference = 5,
        type = "all"
      )

      act_diff_data <- getData(diffMethData)
      act_diff_data <- as.data.frame(act_diff_data, stringsAsFactors = FALSE)

      if (nrow(act_diff_data) == 0) {
        warning(
          "No differential methylation rows returned for ",
          experiment,
          " / ",
          context
        )
      }

      act_diff_data$p_fdr <- p.adjust(
        act_diff_data$pvalue,
        method = "fdr"
      )

      act_diff_data_matched <- add_annotation_chr_column(
        diff_df = act_diff_data,
        methylation_matching = methylation_matching
      )

      raw_output_file <- file.path(
        model_output_dir,
        paste0(
          "Actinia_DMBs_d5_",
          experiment,
          "_",
          context,
          "_covariates_no_outlier_pairs_raw.txt"
        )
      )

      write.table(
        act_diff_data_matched,
        file = raw_output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
      )

      act_diff_df_ann <- act_diff_data_matched

      if (!is.null(genome)) {
        act_diff_df_ann <- tryCatch({
          annotate_diff_table(
            diff_df = act_diff_data_matched,
            genome = genome,
            seq_col = "chr_for_annotation"
          )
        }, error = function(e) {
          message(
            "Warning: GFF3 annotation failed for ",
            experiment,
            " / ",
            context,
            ": ",
            e$message
          )
          message("Continuing with unannotated differential methylation output.")
          act_diff_data_matched
        })
      }

      annotated_output_file <- file.path(
        model_output_dir,
        paste0(
          "act_diff_df_ann_",
          experiment,
          "_",
          context,
          "_covariates_no_outlier_pairs.txt"
        )
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
      message(
        "Error in differential methylation processing for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
    })

    if (is.null(act_diff_data) || is.null(act_diff_df_ann)) {
      gc()
      next
    }

    # -------------------------------
    # Volcano plot
    # -------------------------------

    tryCatch({
      volcano_data <- as.data.frame(act_diff_data, stringsAsFactors = FALSE)

      volcano_data <- volcano_data %>%
        mutate(
          meth.diff = as.numeric(meth.diff),
          p_fdr = as.numeric(p_fdr),
          p_fdr_plot = pmax(p_fdr, .Machine$double.xmin),
          neg_log10_fdr = -log10(p_fdr_plot),
          diffmeth = case_when(
            meth.diff >= 5 ~ "UP",
            meth.diff <= -5 ~ "DOWN",
            TRUE ~ "Not significant"
          ),
          diffmeth = factor(
            diffmeth,
            levels = c("DOWN", "Not significant", "UP")
          )
        ) %>%
        filter(
          is.finite(meth.diff),
          is.finite(neg_log10_fdr)
        )

      if (nrow(volcano_data) == 0) {
        stop("No finite rows available for volcano plot.")
      }

      volcano_plot <- ggplot(
        data = volcano_data,
        aes(x = meth.diff, y = neg_log10_fdr, col = diffmeth)
      ) +
        geom_point(alpha = 0.7, size = 1) +
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
          values = c(
            "DOWN" = "deepskyblue",
            "Not significant" = "grey70",
            "UP" = "brown1"
          ),
          labels = c(
            "Down methylated",
            "Not significant",
            "Up methylated"
          )
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
            "covariate model, outlier pairs excluded"
          )
        ) +
        coord_cartesian(
          xlim = c(-40, 40),
          ylim = c(0, 3)
        ) +
        scale_x_continuous(
          breaks = seq(-40, 40, by = 5)
        ) +
        scale_y_continuous(
          breaks = seq(0, 3, by = 0.2)
        )

      filename <- file.path(
        volcano_figure_dir,
        paste0(
          "volcano_plot_",
          experiment,
          "_",
          context,
          "_covariates_no_outlier_pairs.png"
        )
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
      message(
        "Error creating volcano plot for ",
        experiment,
        " / ",
        context,
        ": ",
        e$message
      )
    })

    rm(
      act_data,
      act_data_f,
      act_data_f_norm,
      act_data_normed_fu,
      dmb_data_exp,
      act_diff_data,
      act_diff_data_matched,
      act_diff_df_ann
    )

    gc()

    message(
      "Completed base-level DMA with covariates: ",
      experiment,
      " / ",
      context
    )
  }
}

message("All base-level cytosine DMA analyses with covariates complete.")