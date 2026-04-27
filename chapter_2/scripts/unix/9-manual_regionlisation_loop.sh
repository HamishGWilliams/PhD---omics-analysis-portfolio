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
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

suppressPackageStartupMessages({
  library(plyranges)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(data.table)
})

# -------------------------------
# Project paths
# -------------------------------

project_dir <- "/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2"

bedgraph_dir <- file.path(project_dir, "data", "processed", "bedgraph")
external_dir <- file.path(project_dir, "data", "external")

# Preferred filename from the RepeatModeler/RepeatMasker + flank workflow
annotation_gff3 <- file.path(external_dir, "combined_annotations.gff3")

# Fallback in case the file was named without the final "s"
if (!file.exists(annotation_gff3)) {
  annotation_gff3_alt <- file.path(external_dir, "combined_annotation.gff3")

  if (file.exists(annotation_gff3_alt)) {
    annotation_gff3 <- annotation_gff3_alt
  }
}

export_dir <- file.path(project_dir, "results", "regionalised_methylation")
overlap_dir <- file.path(export_dir, "overlap_info")
count_dir <- file.path(export_dir, "aggregated_counts")
log_dir <- file.path(project_dir, "logs")

dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(overlap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(count_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(export_dir, "process_log.txt")

cat("Regionalised methylation process log\n", file = log_file)
cat("Started:", as.character(Sys.time()), "\n\n", file = log_file, append = TRUE)

# -------------------------------
# Helper functions
# -------------------------------

log_message <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

get_mcol <- function(gr, col, default = NA_character_) {
  meta <- as.data.frame(mcols(gr))

  if (col %in% names(meta)) {
    as.character(meta[[col]])
  } else {
    rep(default, length(gr))
  }
}

read_bismark_cov_as_granges <- function(file_path) {
  # Bismark coverage files are expected to contain:
  # chr, start, end, methylation_percentage, methylated_count, unmethylated_count

  cov_dt <- fread(
    file_path,
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "percent_methylation", "numCs", "numTs")
  )

  if (nrow(cov_dt) == 0) {
    stop("Coverage file is empty: ", file_path)
  }

  cov_dt[, coverage := numCs + numTs]

  cov_dt <- cov_dt[
    !is.na(chr) &
      !is.na(start) &
      !is.na(end) &
      !is.na(coverage) &
      !is.na(numCs) &
      !is.na(numTs)
  ]

  if (nrow(cov_dt) == 0) {
    stop("No valid rows found after filtering: ", file_path)
  }

  GRanges(
    seqnames = cov_dt$chr,
    ranges = IRanges(start = cov_dt$start, end = cov_dt$end),
    strand = "*",
    coverage = cov_dt$coverage,
    numCs = cov_dt$numCs,
    numTs = cov_dt$numTs,
    percent_methylation = cov_dt$percent_methylation
  )
}

make_empty_region_table <- function(genome_type) {
  data.table(
    feature_id = get_mcol(genome_type, "ID"),
    parent_id = get_mcol(genome_type, "Parent"),
    feature_name = get_mcol(genome_type, "Name"),
    feature_type = get_mcol(genome_type, "type"),
    chr = as.character(seqnames(genome_type)),
    start = start(genome_type),
    end = end(genome_type),
    strand = as.character(strand(genome_type)),
    coverage = 0,
    numCs = 0,
    numTs = 0
  )
}

aggregate_methylation_by_regions <- function(methyl_gr, genome_type) {
  common_seqlevels <- intersect(seqlevels(methyl_gr), seqlevels(genome_type))

  if (length(common_seqlevels) == 0) {
    return(list(
      results_dt = make_empty_region_table(genome_type),
      overlaps = NULL,
      n_overlaps = 0,
      genome_type = genome_type,
      methyl_gr = methyl_gr
    ))
  }

  methyl_gr <- keepSeqlevels(
    methyl_gr,
    common_seqlevels,
    pruning.mode = "coarse"
  )

  genome_type <- keepSeqlevels(
    genome_type,
    common_seqlevels,
    pruning.mode = "coarse"
  )

  results_dt <- make_empty_region_table(genome_type)

  overlaps <- findOverlaps(
    query = genome_type,
    subject = methyl_gr,
    ignore.strand = TRUE
  )

  n_overlaps <- length(overlaps)

  if (n_overlaps > 0) {
    region_idx <- queryHits(overlaps)
    methyl_idx <- subjectHits(overlaps)

    overlap_dt <- data.table(
      region_id = region_idx,
      coverage = mcols(methyl_gr)$coverage[methyl_idx],
      numCs = mcols(methyl_gr)$numCs[methyl_idx],
      numTs = mcols(methyl_gr)$numTs[methyl_idx]
    )

    agg_dt <- overlap_dt[
      ,
      .(
        coverage = sum(coverage, na.rm = TRUE),
        numCs = sum(numCs, na.rm = TRUE),
        numTs = sum(numTs, na.rm = TRUE)
      ),
      by = region_id
    ]

    results_dt[agg_dt$region_id, coverage := agg_dt$coverage]
    results_dt[agg_dt$region_id, numCs := agg_dt$numCs]
    results_dt[agg_dt$region_id, numTs := agg_dt$numTs]
  }

  list(
    results_dt = results_dt,
    overlaps = overlaps,
    n_overlaps = n_overlaps,
    genome_type = genome_type,
    methyl_gr = methyl_gr
  )
}

# -------------------------------
# Load combined annotation file
# -------------------------------

if (!file.exists(annotation_gff3)) {
  stop(
    "Combined annotation GFF3 not found. Expected one of:\n",
    file.path(external_dir, "combined_annotations.gff3"), "\n",
    file.path(external_dir, "combined_annotation.gff3")
  )
}

log_message("Loading combined annotation file: ", annotation_gff3)

genome <- read_gff3(annotation_gff3)

# Convert through data.frame to make metadata columns predictable and remove NA values
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

if (!"type" %in% names(genome_df)) {
  stop("The combined GFF3 does not contain a usable 'type' column.")
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

types <- unique(as.character(mcols(genome)$type))
types <- types[!is.na(types) & types != "" & types != "NA"]

# Put the biologically relevant methylation-region classes first if present.
preferred_types <- c(
  "gene",
  "exon",
  "five_prime_UTR",
  "three_prime_UTR",
  "promoter",
  "downstream_region",
  "dispersed_repeat",
  "ncRNA_gene"
)

types <- c(
  intersect(preferred_types, types),
  setdiff(types, preferred_types)
)

log_message("Detected ", length(types), " annotation feature types:")
log_message(paste(types, collapse = ", "))

if (!"promoter" %in% types) {
  log_message("WARNING: No promoter features detected in combined annotation file.")
}

if (!"downstream_region" %in% types) {
  log_message("WARNING: No downstream_region features detected in combined annotation file.")
}

if (!"dispersed_repeat" %in% types) {
  log_message("WARNING: No dispersed_repeat features detected in combined annotation file.")
}

# -------------------------------
# Sample definitions
# -------------------------------

contexts <- c("CpG", "CHG", "CHH")

experiments <- list(
  exp1 = data.frame(
    sample_dir = c(
      "Sample_3-3_D",
      "Sample_9-9_D",
      "Sample_13-13_D",
      "Sample_15-15_D",
      "Sample_16-16_D",
      "Sample_17-17_D",
      "Sample_22-22_D",
      "Sample_27-27_D"
    ),
    sample_id = c("3", "9", "13", "15", "16", "17", "22", "27"),
    treatment = c(1, 0, 0, 1, 0, 0, 1, 1),
    stringsAsFactors = FALSE
  ),
  exp2 = data.frame(
    sample_dir = c(
      "Sample_6-6_D",
      "Sample_7-7_D",
      "Sample_8-8_D",
      "Sample_10-10_D",
      "Sample_14-14_D",
      "Sample_20-20_D",
      "Sample_21-21_D",
      "Sample_29-11_redo_D",
      "Sample_32-18_redo_D",
      "Sample_36-2_redo_D"
    ),
    sample_id = c("6", "7", "8", "10", "14", "20", "21", "29_11", "32_18", "36_2"),
    treatment = c(0, 0, 1, 0, 1, 1, 1, 0, 0, 1),
    stringsAsFactors = FALSE
  )
)

# -------------------------------
# Main regionalisation loop
# -------------------------------

for (exp_name in names(experiments)) {
  sample_info <- experiments[[exp_name]]

  log_message("Starting experiment: ", exp_name)

  for (context in contexts) {
    log_message("  Starting context: ", context)

    for (sample_idx in seq_len(nrow(sample_info))) {
      sample_dir <- sample_info$sample_dir[sample_idx]
      sample_id <- sample_info$sample_id[sample_idx]
      treatment <- sample_info$treatment[sample_idx]

      cov_file <- file.path(
        bedgraph_dir,
        sample_dir,
        paste0("meth_", context, "_cov_reads")
      )

      if (!file.exists(cov_file)) {
        log_message(
          "    Missing file for ",
          exp_name,
          " / ",
          context,
          " / ",
          sample_id,
          ": ",
          cov_file
        )
        next
      }

      log_message("    Reading sample ", sample_id, " from ", cov_file)

      methyl_gr <- tryCatch(
        read_bismark_cov_as_granges(cov_file),
        error = function(e) {
          log_message("    Error reading ", cov_file, ": ", e$message)
          NULL
        }
      )

      if (is.null(methyl_gr)) {
        next
      }

      for (current_type in types) {
        safe_type <- safe_name(current_type)

        log_message("      Processing feature type: ", current_type)

        genome_type <- genome[mcols(genome)$type == current_type]

        if (length(genome_type) == 0) {
          log_message("      No regions found for feature type: ", current_type)
          next
        }

        regionalised <- tryCatch(
          aggregate_methylation_by_regions(
            methyl_gr = methyl_gr,
            genome_type = genome_type
          ),
          error = function(e) {
            log_message(
              "      Error aggregating sample ",
              sample_id,
              " / ",
              context,
              " / ",
              current_type,
              ": ",
              e$message
            )
            NULL
          }
        )

        if (is.null(regionalised)) {
          next
        }

        log_message(
          "      Number of overlaps for sample ",
          sample_id,
          " / ",
          context,
          " / ",
          current_type,
          ": ",
          regionalised$n_overlaps
        )

        # -------------------------------
        # Save overlap information
        # -------------------------------

        overlap_output_file <- file.path(
          overlap_dir,
          exp_name,
          context,
          safe_type,
          paste0(sample_id, "_", safe_type, "_overlaps_info.tsv")
        )

        dir.create(dirname(overlap_output_file), recursive = TRUE, showWarnings = FALSE)

        if (!is.null(regionalised$overlaps) && regionalised$n_overlaps > 0) {
          overlaps <- regionalised$overlaps
          genome_type_filtered <- regionalised$genome_type
          methyl_gr_filtered <- regionalised$methyl_gr

          query_gr <- genome_type_filtered[queryHits(overlaps)]
          subject_gr <- methyl_gr_filtered[subjectHits(overlaps)]

          overlaps_info <- data.frame(
            experiment = exp_name,
            context = context,
            sample_id = sample_id,
            sample_dir = sample_dir,
            treatment = treatment,
            feature_type = current_type,
            feature_id = get_mcol(query_gr, "ID"),
            parent_id = get_mcol(query_gr, "Parent"),
            feature_name = get_mcol(query_gr, "Name"),
            region_chr = as.character(seqnames(query_gr)),
            region_start = start(query_gr),
            region_end = end(query_gr),
            region_strand = as.character(strand(query_gr)),
            methylation_chr = as.character(seqnames(subject_gr)),
            methylation_position = start(subject_gr),
            coverage = mcols(subject_gr)$coverage,
            numCs = mcols(subject_gr)$numCs,
            numTs = mcols(subject_gr)$numTs,
            percent_methylation = mcols(subject_gr)$percent_methylation
          )
        } else {
          overlaps_info <- data.frame(
            experiment = character(),
            context = character(),
            sample_id = character(),
            sample_dir = character(),
            treatment = numeric(),
            feature_type = character(),
            feature_id = character(),
            parent_id = character(),
            feature_name = character(),
            region_chr = character(),
            region_start = integer(),
            region_end = integer(),
            region_strand = character(),
            methylation_chr = character(),
            methylation_position = integer(),
            coverage = numeric(),
            numCs = numeric(),
            numTs = numeric(),
            percent_methylation = numeric()
          )
        }

        write.table(
          overlaps_info,
          file = overlap_output_file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )

        # -------------------------------
        # Save aggregated regional counts
        # -------------------------------

        count_output_file <- file.path(
          count_dir,
          exp_name,
          context,
          safe_type,
          paste0(sample_id, "_", safe_type, "_aggregated_methylation_counts.txt")
        )

        dir.create(dirname(count_output_file), recursive = TRUE, showWarnings = FALSE)

        results_dt <- copy(regionalised$results_dt)

        results_dt[, experiment := exp_name]
        results_dt[, context := context]
        results_dt[, sample_id := sample_id]
        results_dt[, sample_dir := sample_dir]
        results_dt[, treatment := treatment]

        results_dt[, percent_methylated := fifelse(
          coverage > 0,
          (numCs / coverage) * 100,
          NA_real_
        )]

        setcolorder(
          results_dt,
          c(
            "experiment",
            "context",
            "sample_id",
            "sample_dir",
            "treatment",
            "feature_type",
            "feature_id",
            "parent_id",
            "feature_name",
            "chr",
            "start",
            "end",
            "strand",
            "coverage",
            "numCs",
            "numTs",
            "percent_methylated"
          )
        )

        fwrite(
          results_dt,
          file = count_output_file,
          sep = "\t",
          quote = FALSE,
          na = "NA"
        )

        log_message("      Saved aggregated counts: ", count_output_file)
      }

      rm(methyl_gr)
      gc()
    }
  }
}

cat("\nProcessing completed: ", as.character(Sys.time()), "\n", file = log_file, append = TRUE)
log_message("Regionalised methylation processing complete.")