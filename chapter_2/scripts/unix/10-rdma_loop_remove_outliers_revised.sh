#!/bin/bash
#SBATCH --job-name=region_level_dma_no_outliers
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

module load r/4.2.2

Rscript --vanilla - <<'RSCRIPT'

suppressPackageStartupMessages({
  library(methylKit)
  library(plyranges)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(grid)
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

metadata_file <- file.path(
  project_dir,
  "data",
  "metadata",
  "Methyl_sample_groups.txt"
)

results_dir <- file.path(project_dir, "results")

model_output_dir <- file.path(
  results_dir,
  "model_outputs",
  "region_level_dma_no_outlier_pairs"
)

figure_dir <- file.path(
  results_dir,
  "figures",
  "region_level_dma_no_outlier_pairs"
)

table_output_dir <- file.path(
  results_dir,
  "tables",
  "region_level_dma_no_outlier_pairs"
)

qc_figure_dir <- file.path(figure_dir, "qc")
pca_figure_dir <- file.path(figure_dir, "pca")
volcano_figure_dir <- file.path(figure_dir, "volcano")
venn_figure_dir <- file.path(figure_dir, "venn")

dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(venn_figure_dir, recursive = TRUE, showWarnings = FALSE)

# Temporary methylKit input files are written to job-local scratch only.
# They are deleted immediately after methRead() has loaded them into memory.

job_tmp_base <- Sys.getenv("SLURM_TMPDIR", unset = tempdir())

tmp_methylkit_dir <- file.path(
  job_tmp_base,
  paste0("region_level_dma_methylkit_inputs_", Sys.getpid())
)

dir.create(tmp_methylkit_dir, recursive = TRUE, showWarnings = FALSE)

on.exit({
  if (dir.exists(tmp_methylkit_dir)) {
    unlink(tmp_methylkit_dir, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

setwd(job_tmp_base)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "24"))
message("Using ", n_cores, " cores.")
message("Temporary methylKit files will be written to: ", tmp_methylkit_dir)

# -------------------------------
# Analysis settings
# -------------------------------

experiments <- c("exp1", "exp2")
contexts <- c("CpG", "CHG", "CHH")

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

fdr_threshold <- 0.1
meth_diff_threshold <- 5

make_qc_plots <- TRUE
make_pca_plots <- TRUE
make_venn_outputs <- TRUE

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

# -------------------------------
# Final output files only
# -------------------------------

combined_rdata_file <- file.path(
  model_output_dir,
  "region_level_dma_no_outlier_pairs_combined_results.RData"
)

combined_txt_file <- file.path(
  model_output_dir,
  "region_level_dma_no_outlier_pairs_combined_results.txt"
)

dmr_count_file <- file.path(
  table_output_dir,
  "region_level_dma_no_outlier_pairs_hyper_hypo_DMR_counts.txt"
)

gene_intersection_file <- file.path(
  table_output_dir,
  "region_level_dma_no_outlier_pairs_associated_gene_context_intersections.txt"
)

dmr_intersection_file <- file.path(
  table_output_dir,
  "region_level_dma_no_outlier_pairs_DMR_context_intersections.txt"
)

volcano_file <- file.path(
  volcano_figure_dir,
  "region_level_dma_no_outlier_pairs_combined_volcano_facet_by_region.png"
)

qc_methylation_pdf <- file.path(
  qc_figure_dir,
  "region_level_dma_no_outlier_pairs_methylation_QC_all.pdf"
)

qc_coverage_pdf <- file.path(
  qc_figure_dir,
  "region_level_dma_no_outlier_pairs_coverage_QC_all.pdf"
)

pca_pdf <- file.path(
  pca_figure_dir,
  "region_level_dma_no_outlier_pairs_PCA_all.pdf"
)

gene_venn_pdf <- file.path(
  venn_figure_dir,
  "region_level_dma_no_outlier_pairs_associated_gene_context_venn.pdf"
)

dmr_venn_pdf <- file.path(
  venn_figure_dir,
  "region_level_dma_no_outlier_pairs_DMR_context_venn.pdf"
)

# -------------------------------
# Plot device helpers
# -------------------------------

plot_header_page <- function(title, subtitle = NULL) {
  plot.new()
  text(0.5, 0.6, title, cex = 1.2, font = 2)

  if (!is.null(subtitle)) {
    text(0.5, 0.45, subtitle, cex = 0.9)
  }
}

qc_methylation_dev <- NULL
qc_coverage_dev <- NULL
pca_dev <- NULL

if (make_qc_plots) {
  pdf(qc_methylation_pdf, width = 10, height = 8)
  qc_methylation_dev <- dev.cur()

  pdf(qc_coverage_pdf, width = 12, height = 8)
  qc_coverage_dev <- dev.cur()
}

if (make_pca_plots) {
  pdf(pca_pdf, width = 9, height = 7)
  pca_dev <- dev.cur()
}

close_open_plot_devices <- function() {
  open_devs <- dev.list()

  if (!is.null(open_devs)) {
    open_dev_ids <- as.integer(open_devs)

    for (dev_id in c(qc_methylation_dev, qc_coverage_dev, pca_dev)) {
      if (!is.null(dev_id) && dev_id %in% open_dev_ids) {
        dev.set(dev_id)
        dev.off()
      }
    }
  }
}

on.exit(close_open_plot_devices(), add = TRUE)

# -------------------------------
# Metadata
# -------------------------------

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

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

genome <- tryCatch(
  {
    read_gff3(annotation_gff3)
  },
  error = function(e) {
    stop(
      "Failed to read GFF3 annotation file: ",
      annotation_gff3,
      "\nError message: ",
      e$message
    )
  }
)

genome_df <- as.data.frame(genome)

genome_df[] <- lapply(genome_df, function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x[is.na(x)] <- "NA"
  }

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

if (is.null(mcols(genome)$type)) {
  stop("The GFF3 annotation does not contain a usable 'type' column.")
}

detected_regions <- unique(as.character(mcols(genome)$type))
detected_regions <- detected_regions[
  !is.na(detected_regions) &
    detected_regions != "" &
    detected_regions != "NA"
]

regions <- intersect(preferred_regions, detected_regions)

if (length(regions) == 0) {
  stop(
    "None of the preferred region types were found in the GFF3. Detected types: ",
    paste(detected_regions, collapse = ", ")
  )
}

message("Region types selected for analysis:")
message(paste(regions, collapse = ", "))

# -------------------------------
# Sample definitions: outlier pairs excluded
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
# Temporary input preparation
# -------------------------------

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
# Annotation helpers
# -------------------------------

clean_gene_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^gene:", "", x)
  x <- gsub("^transcript:", "", x)
  x <- gsub("^ID=", "", x)
  x <- gsub("^Parent=", "", x)
  x[x %in% c("", "NA", ".", "nan", "NaN")] <- NA_character_
  x
}

get_annotation_gene_names <- function(gr) {
  m <- as.data.frame(mcols(gr), stringsAsFactors = FALSE)

  candidate_cols <- c(
    "gene_name",
    "gene",
    "Name",
    "name",
    "gene_id",
    "locus_tag",
    "ID",
    "Parent"
  )

  existing_cols <- intersect(candidate_cols, names(m))

  if (length(existing_cols) == 0) {
    return(rep(NA_character_, length(gr)))
  }

  gene_name <- rep(NA_character_, nrow(m))

  for (col in existing_cols) {
    values <- clean_gene_name(m[[col]])
    replace_idx <- is.na(gene_name) | gene_name == ""
    gene_name[replace_idx] <- values[replace_idx]
  }

  clean_gene_name(gene_name)
}

annotate_diff_table <- function(diff_df, genome_region, experiment, context, region) {
  out <- as.data.table(diff_df)

  out[, experiment := experiment]
  out[, context := context]
  out[, region_type := region]

  out[, chr := as.character(chr)]
  out[, start := as.integer(start)]
  out[, end := as.integer(end)]
  out[, meth.diff := as.numeric(meth.diff)]
  out[, pvalue := as.numeric(pvalue)]

  out[, p_fdr := p.adjust(pvalue, method = "fdr")]

  out[, dmr_status := fifelse(
    p_fdr <= fdr_threshold & meth.diff >= meth_diff_threshold,
    "Hypermethylated",
    fifelse(
      p_fdr <= fdr_threshold & meth.diff <= -meth_diff_threshold,
      "Hypomethylated",
      "Not significant"
    )
  )]

  out[, result_row_id := seq_len(.N)]
  out[, annotation_gene_name := NA_character_]
  out[, annotation_overlap_count := 0L]

  valid <- !is.na(out$chr) &
    !is.na(out$start) &
    !is.na(out$end) &
    out$start <= out$end

  if (nrow(out) > 0 && any(valid) && length(genome_region) > 0) {
    valid_idx <- which(valid)

    query_gr <- GRanges(
      seqnames = out$chr[valid],
      ranges = IRanges(
        start = out$start[valid],
        end = out$end[valid]
      ),
      strand = "*"
    )

    hits <- suppressWarnings(
      findOverlaps(
        query_gr,
        genome_region,
        ignore.strand = TRUE
      )
    )

    if (length(hits) > 0) {
      hit_query_rows <- valid_idx[queryHits(hits)]
      hit_subject_rows <- subjectHits(hits)

      subject_gene_names <- get_annotation_gene_names(genome_region)

      hit_dt <- data.table(
        result_row_id = hit_query_rows,
        annotation_gene_name = subject_gene_names[hit_subject_rows]
      )

      overlap_count <- hit_dt[, .(
        annotation_overlap_count_new = .N
      ), by = result_row_id]

      gene_summary <- hit_dt[
        !is.na(annotation_gene_name) &
          annotation_gene_name != "",
        .(
          annotation_gene_name_new = paste(
            sort(unique(annotation_gene_name)),
            collapse = ";"
          )
        ),
        by = result_row_id
      ]

      out <- merge(
        out,
        overlap_count,
        by = "result_row_id",
        all.x = TRUE,
        sort = FALSE
      )

      out <- merge(
        out,
        gene_summary,
        by = "result_row_id",
        all.x = TRUE,
        sort = FALSE
      )

      out[
        !is.na(annotation_overlap_count_new),
        annotation_overlap_count := annotation_overlap_count_new
      ]

      out[
        !is.na(annotation_gene_name_new),
        annotation_gene_name := annotation_gene_name_new
      ]

      out[, annotation_overlap_count_new := NULL]
      out[, annotation_gene_name_new := NULL]
    }
  }

  out[, dmr_id := paste(region_type, chr, start, end, sep = "|")]

  out[, dmr_gene_region_id := paste(
    region_type,
    ifelse(
      !is.na(annotation_gene_name) & annotation_gene_name != "",
      annotation_gene_name,
      "unannotated"
    ),
    chr,
    start,
    end,
    sep = "|"
  )]

  setorder(out, result_row_id)
  out[, result_row_id := NULL]

  out[]
}

# -------------------------------
# PCA helper
# -------------------------------

plot_regional_pca <- function(
  act_data_normed_fu,
  experiment,
  context,
  region,
  sample_metadata
) {
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
    slice(chull(PC1, PC2)) %>%
    ungroup()

  ggplot(pca_data, aes(x = PC1, y = PC2)) +
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
        "Regional PCA:",
        context,
        region,
        experiment,
        "outlier pairs excluded"
      ),
      x = paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)"),
      color = "Group",
      fill = "Group"
    )
}

# -------------------------------
# Main analysis loop
# -------------------------------

all_results_list <- list()

for (experiment in experiments) {
  sample_info <- get_sample_info(experiment)

  for (context in contexts) {
    for (region in regions) {
      safe_region <- safe_name(region)

      message("------------------------------------------------------------")
      message("Starting region-level DMA")
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
        unlink(paths$tmp_files, force = TRUE)
        next
      }

      file.list <- as.list(converted_files)
      sample.id <- as.list(sample_info$methylkit_sample_id)
      treatment <- sample_info$treatment
      covariates_df <- data.frame(ind_id = factor(sample_info$ind_id))

      # -------------------------------
      # Read methylation data into memory
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

      # Temporary input files are no longer needed after methRead().
      unlink(paths$tmp_files, force = TRUE)

      if (is.null(act_data)) {
        gc()
        next
      }

      # -------------------------------
      # QC plots: combined PDFs only
      # -------------------------------

      if (make_qc_plots) {
        tryCatch({
          dev.set(qc_methylation_dev)
          plot_header_page(
            paste("Methylation QC:", experiment, context, region),
            "Each following page is one sample."
          )
          lapply(act_data, getMethylationStats, plot = TRUE)
        }, error = function(e) {
          message("Error plotting methylation stats: ", e$message)
        })

        tryCatch({
          dev.set(qc_coverage_dev)
          plot_header_page(
            paste("Coverage QC:", experiment, context, region),
            "Each following page is one sample."
          )
          lapply(act_data, getCoverageStats, plot = TRUE)
        }, error = function(e) {
          message("Error plotting coverage stats: ", e$message)
        })
      }

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
          suffix = paste0(context, "_", safe_region, "_f")
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
          suffix = paste0(context, "_", safe_region, "_normed_fu3")
        )
      }, error = function(e) {
        message("Error uniting methylation data: ", e$message)
      })

      if (is.null(act_data_normed_fu)) {
        gc()
        next
      }

      # -------------------------------
      # PCA: combined PDF only
      # -------------------------------

      if (make_pca_plots) {
        tryCatch({
          pca_plot <- plot_regional_pca(
            act_data_normed_fu = act_data_normed_fu,
            experiment = experiment,
            context = context,
            region = region,
            sample_metadata = sample_metadata
          )

          dev.set(pca_dev)
          print(pca_plot)
        }, error = function(e) {
          message("Error plotting PCA: ", e$message)
        })
      }

      # -------------------------------
      # Differential methylation
      # -------------------------------

      dmb_data_exp <- NULL

      tryCatch({
        dmb_data_exp <- calculateDiffMeth(
          act_data_normed_fu,
          covariates = covariates_df,
          overdispersion = "MN",
          mc.cores = n_cores,
          suffix = paste0(
            context,
            "_fu3_",
            safe_region,
            "_odMN_region_level_no_outlier_pairs"
          )
        )
      }, error = function(e) {
        message("Error in differential methylation calculation: ", e$message)
      })

      if (is.null(dmb_data_exp)) {
        gc()
        next
      }

      # -------------------------------
      # Extract, add FDR, annotate
      # -------------------------------

      tryCatch({
        act_diff_data <- getData(dmb_data_exp)
        act_diff_data <- as.data.frame(act_diff_data, stringsAsFactors = FALSE)

        genome_region <- genome[mcols(genome)$type == region]

        act_diff_data_ann <- annotate_diff_table(
          diff_df = act_diff_data,
          genome_region = genome_region,
          experiment = experiment,
          context = context,
          region = region
        )

        result_key <- paste(experiment, context, safe_region, sep = "__")
        all_results_list[[result_key]] <- act_diff_data_ann

        message(
          "Stored results in memory: ",
          experiment,
          " / ",
          context,
          " / ",
          region,
          " rows = ",
          nrow(act_diff_data_ann)
        )

      }, error = function(e) {
        message(
          "Error extracting/annotating differential methylation results for ",
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

# Close multipage plot PDFs before final outputs
close_open_plot_devices()

# -------------------------------
# Combine results
# -------------------------------

if (length(all_results_list) == 0) {
  stop("No differential methylation results were generated.")
}

all_results <- rbindlist(
  all_results_list,
  use.names = TRUE,
  fill = TRUE
)

all_results[, experiment := factor(experiment, levels = experiments)]
all_results[, context := factor(context, levels = contexts)]
all_results[, region_type := factor(region_type, levels = regions)]

message("Combined result rows: ", nrow(all_results))

# -------------------------------
# Save combined results
# -------------------------------

fwrite(
  all_results,
  file = combined_txt_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

message("Saved combined TXT results: ", combined_txt_file)

# -------------------------------
# Hyper/hypo DMR count table
# -------------------------------

dmr_count_table <- all_results[, .(
  total_tested = .N,
  hyper_methylated = sum(dmr_status == "Hypermethylated", na.rm = TRUE),
  hypo_methylated = sum(dmr_status == "Hypomethylated", na.rm = TRUE),
  total_significant_DMRs = sum(
    dmr_status %in% c("Hypermethylated", "Hypomethylated"),
    na.rm = TRUE
  )
), by = .(
  experiment,
  region_type,
  context
)]

setorder(dmr_count_table, experiment, region_type, context)

fwrite(
  dmr_count_table,
  file = dmr_count_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

message("Saved DMR count table: ", dmr_count_file)

# -------------------------------
# Single combined volcano plot
# -------------------------------

volcano_data <- copy(all_results)

volcano_data[, pvalue := as.numeric(pvalue)]
volcano_data[, meth.diff := as.numeric(meth.diff)]
volcano_data[, pvalue_plot := pmax(pvalue, .Machine$double.xmin)]
volcano_data[, neg_log10_pvalue := -log10(pvalue_plot)]

volcano_data <- volcano_data[
  is.finite(meth.diff) &
    is.finite(neg_log10_pvalue)
]

volcano_data[, dmr_status := factor(
  dmr_status,
  levels = c("Hypomethylated", "Not significant", "Hypermethylated")
)]

combined_volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = meth.diff,
    y = neg_log10_pvalue,
    colour = dmr_status,
    shape = context
  )
) +
  geom_point(alpha = 0.65, size = 0.7) +
  geom_vline(
    xintercept = -meth_diff_threshold,
    colour = "blue",
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = meth_diff_threshold,
    colour = "red",
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(fdr_threshold),
    colour = "black",
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ region_type,
    ncol = 2
  ) +
  scale_colour_manual(
    values = c(
      "Hypomethylated" = "deepskyblue",
      "Not significant" = "grey70",
      "Hypermethylated" = "brown1"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "CpG" = 16,
      "CHG" = 17,
      "CHH" = 15
    ),
    drop = FALSE
  ) +
  coord_cartesian(
    xlim = c(-40, 40),
    ylim = c(0, 6)
  ) +
  scale_x_continuous(
    breaks = seq(-40, 40, by = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 1)
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(
      face = "bold",
      margin = margin(0, 20, 0, 0),
      colour = "black"
    ),
    axis.title.x = element_text(
      hjust = 0.5,
      face = "bold",
      margin = margin(20, 0, 0, 0),
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = rel(1.05)
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "Differential methylation %",
    y = expression("-log"[10] * " raw p-value"),
    colour = "DMR status",
    shape = "Cytosine context"
  ) +
  ggtitle(
    "Region-level differential methylation across annotation regions"
  )

ggsave(
  filename = volcano_file,
  plot = combined_volcano_plot,
  width = 24,
  height = 32,
  units = "cm",
  dpi = 300
)

message("Saved combined volcano plot: ", volcano_file)

# -------------------------------
# Intersection and Venn helpers
# -------------------------------

make_gene_long_table <- function(dt) {
  sig_dt <- dt[
    dmr_status %in% c("Hypermethylated", "Hypomethylated") &
      !is.na(annotation_gene_name) &
      annotation_gene_name != ""
  ]

  if (nrow(sig_dt) == 0) {
    return(data.table(
      experiment = character(),
      region_type = character(),
      context = character(),
      entity_type = character(),
      entity_id = character(),
      entity_label = character()
    ))
  }

  rows <- lapply(seq_len(nrow(sig_dt)), function(i) {
    genes <- unlist(
      strsplit(
        as.character(sig_dt$annotation_gene_name[i]),
        ";",
        fixed = TRUE
      )
    )

    genes <- unique(trimws(genes))
    genes <- genes[
      !is.na(genes) &
        genes != "" &
        genes != "NA"
    ]

    if (length(genes) == 0) {
      return(NULL)
    }

    data.table(
      experiment = as.character(sig_dt$experiment[i]),
      region_type = as.character(sig_dt$region_type[i]),
      context = as.character(sig_dt$context[i]),
      entity_type = "associated_gene",
      entity_id = genes,
      entity_label = genes
    )
  })

  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

make_dmr_long_table <- function(dt) {
  sig_dt <- dt[
    dmr_status %in% c("Hypermethylated", "Hypomethylated")
  ]

  if (nrow(sig_dt) == 0) {
    return(data.table(
      experiment = character(),
      region_type = character(),
      context = character(),
      entity_type = character(),
      entity_id = character(),
      entity_label = character()
    ))
  }

  unique(sig_dt[, .(
    experiment = as.character(experiment),
    region_type = as.character(region_type),
    context = as.character(context),
    entity_type = "DMR",
    entity_id = dmr_id,
    entity_label = dmr_gene_region_id
  )])
}

make_intersection_table <- function(entity_long, entity_type_label) {
  if (nrow(entity_long) == 0) {
    return(data.table(
      experiment = character(),
      region_type = character(),
      entity_type = character(),
      entity_id = character(),
      entity_label = character(),
      CpG = integer(),
      CHG = integer(),
      CHH = integer(),
      intersection_group = character()
    ))
  }

  entity_long <- unique(entity_long)

  labels <- entity_long[, .(
    entity_label = paste(sort(unique(entity_label)), collapse = ";")
  ), by = .(
    experiment,
    region_type,
    entity_type,
    entity_id
  )]

  present <- unique(entity_long[, .(
    experiment,
    region_type,
    entity_type,
    entity_id,
    context
  )])

  present[, present := 1L]

  wide <- dcast(
    present,
    experiment + region_type + entity_type + entity_id ~ context,
    value.var = "present",
    fun.aggregate = sum,
    fill = 0
  )

  for (ctx in contexts) {
    if (!ctx %in% names(wide)) {
      wide[, (ctx) := 0L]
    }
  }

  wide[, (contexts) := lapply(.SD, function(x) as.integer(x > 0)),
       .SDcols = contexts]

  wide <- merge(
    wide,
    labels,
    by = c("experiment", "region_type", "entity_type", "entity_id"),
    all.x = TRUE,
    sort = FALSE
  )

  wide[, intersection_group := apply(
    .SD,
    1,
    function(x) {
      active <- contexts[as.logical(as.integer(x))]

      if (length(active) == 0) {
        "none"
      } else {
        paste(active, collapse = "&")
      }
    }
  ), .SDcols = contexts]

  wide[, entity_type := entity_type_label]

  setcolorder(
    wide,
    c(
      "experiment",
      "region_type",
      "entity_type",
      "entity_id",
      "entity_label",
      contexts,
      "intersection_group"
    )
  )

  setorder(wide, experiment, region_type, intersection_group, entity_id)

  wide[]
}

get_venn_counts <- function(set_list) {
  A <- set_list[[contexts[1]]]
  B <- set_list[[contexts[2]]]
  C <- set_list[[contexts[3]]]

  list(
    A_only = length(setdiff(A, union(B, C))),
    B_only = length(setdiff(B, union(A, C))),
    C_only = length(setdiff(C, union(A, B))),
    AB_only = length(setdiff(intersect(A, B), C)),
    AC_only = length(setdiff(intersect(A, C), B)),
    BC_only = length(setdiff(intersect(B, C), A)),
    ABC = length(Reduce(intersect, list(A, B, C))),
    total = length(unique(c(A, B, C)))
  )
}

plot_three_set_venn <- function(set_list, title, subtitle) {
  counts <- get_venn_counts(set_list)

  grid.newpage()

  if (counts$total == 0) {
    grid.text(
      title,
      x = 0.5,
      y = 0.65,
      gp = gpar(fontsize = 14, fontface = "bold")
    )

    grid.text(
      subtitle,
      x = 0.5,
      y = 0.58,
      gp = gpar(fontsize = 10)
    )

    grid.text(
      "No significant DMRs available for this comparison.",
      x = 0.5,
      y = 0.45,
      gp = gpar(fontsize = 11)
    )

    return(invisible(NULL))
  }

  grid.circle(
    x = 0.42,
    y = 0.56,
    r = 0.25,
    gp = gpar(
      fill = adjustcolor("red", alpha.f = 0.18),
      col = "red",
      lwd = 2
    )
  )

  grid.circle(
    x = 0.58,
    y = 0.56,
    r = 0.25,
    gp = gpar(
      fill = adjustcolor("blue", alpha.f = 0.18),
      col = "blue",
      lwd = 2
    )
  )

  grid.circle(
    x = 0.50,
    y = 0.39,
    r = 0.25,
    gp = gpar(
      fill = adjustcolor("green3", alpha.f = 0.18),
      col = "green3",
      lwd = 2
    )
  )

  grid.text(
    title,
    x = 0.5,
    y = 0.94,
    gp = gpar(fontsize = 14, fontface = "bold")
  )

  grid.text(
    subtitle,
    x = 0.5,
    y = 0.89,
    gp = gpar(fontsize = 10)
  )

  grid.text(contexts[1], x = 0.29, y = 0.78, gp = gpar(fontsize = 11, fontface = "bold"))
  grid.text(contexts[2], x = 0.71, y = 0.78, gp = gpar(fontsize = 11, fontface = "bold"))
  grid.text(contexts[3], x = 0.50, y = 0.10, gp = gpar(fontsize = 11, fontface = "bold"))

  grid.text(counts$A_only, x = 0.34, y = 0.58, gp = gpar(fontsize = 12))
  grid.text(counts$B_only, x = 0.66, y = 0.58, gp = gpar(fontsize = 12))
  grid.text(counts$C_only, x = 0.50, y = 0.29, gp = gpar(fontsize = 12))

  grid.text(counts$AB_only, x = 0.50, y = 0.64, gp = gpar(fontsize = 12))
  grid.text(counts$AC_only, x = 0.41, y = 0.43, gp = gpar(fontsize = 12))
  grid.text(counts$BC_only, x = 0.59, y = 0.43, gp = gpar(fontsize = 12))
  grid.text(counts$ABC, x = 0.50, y = 0.50, gp = gpar(fontsize = 12, fontface = "bold"))

  grid.text(
    paste("Total unique:", counts$total),
    x = 0.5,
    y = 0.03,
    gp = gpar(fontsize = 10)
  )

  invisible(NULL)
}

write_venn_pdf <- function(entity_long, output_pdf, entity_label) {
  pdf(output_pdf, width = 8, height = 8)

  for (experiment_i in experiments) {
    for (region_i in regions) {
      subset_dt <- entity_long[
        experiment == experiment_i &
          region_type == region_i
      ]

      set_list <- setNames(
        lapply(contexts, function(ctx) {
          unique(subset_dt[context == ctx, entity_id])
        }),
        contexts
      )

      plot_three_set_venn(
        set_list = set_list,
        title = paste(entity_label, "context overlap"),
        subtitle = paste(experiment_i, "-", region_i)
      )
    }
  }

  dev.off()
}

# -------------------------------
# Venn diagrams and intersection tables
# -------------------------------

gene_long <- make_gene_long_table(all_results)
dmr_long <- make_dmr_long_table(all_results)

gene_intersection_table <- make_intersection_table(
  gene_long,
  entity_type_label = "associated_gene"
)

dmr_intersection_table <- make_intersection_table(
  dmr_long,
  entity_type_label = "DMR"
)

fwrite(
  gene_intersection_table,
  file = gene_intersection_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  dmr_intersection_table,
  file = dmr_intersection_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

message("Saved associated gene intersection table: ", gene_intersection_file)
message("Saved DMR intersection table: ", dmr_intersection_file)

if (make_venn_outputs) {
  write_venn_pdf(
    entity_long = gene_long,
    output_pdf = gene_venn_pdf,
    entity_label = "Associated gene"
  )

  write_venn_pdf(
    entity_long = dmr_long,
    output_pdf = dmr_venn_pdf,
    entity_label = "DMR"
  )

  message("Saved associated gene Venn PDF: ", gene_venn_pdf)
  message("Saved DMR Venn PDF: ", dmr_venn_pdf)
}

# -------------------------------
# Save final RData
# -------------------------------

save(
  all_results,
  dmr_count_table,
  gene_intersection_table,
  dmr_intersection_table,
  file = combined_rdata_file
)

message("Saved combined RData: ", combined_rdata_file)

message("All region-level DMA analyses with outlier pairs excluded are complete.")

RSCRIPT