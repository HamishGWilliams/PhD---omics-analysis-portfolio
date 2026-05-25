library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(knitr)

# -----------------------------
# 1. Load count data
# -----------------------------

counts <- as.matrix(read.table("C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/Processed/A_Equina_Counts_unstranded_copy.txt", header = T, sep = '\t'))

countData <- as.matrix(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)))
{apply(countData, 2, as.numeric)
  sapply(countData, as.numeric)
  class(countData) <- "numeric"
  storage.mode(countData) <- "numeric"}

storage.mode(countData) <- "numeric"

# -----------------------------
# 2. Load sample metadata
# -----------------------------

colData <- read.table(
  "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/MS_colData.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = TRUE,
  check.names = FALSE
)

reorder_idx <- match(colnames(countData), rownames(colData))

if (anyNA(reorder_idx)) {
  stop(
    "Name mismatch after match(); check for typos/duplicates in sample names.\n",
    "Problem samples: ",
    paste(colnames(countData)[is.na(reorder_idx)], collapse = ", ")
  )
}

colData <- colData[reorder_idx, , drop = FALSE]

colData$Diesel   <- factor(colData$Diesel,   levels = c("N", "Y"))
colData$Salinity <- factor(colData$Salinity, levels = c("N", "Y"))

# -----------------------------
# 3. Count-derived QC metrics
# -----------------------------
# These are available directly from the featureCounts count table.
# They represent assigned reads/fragments, not total sequenced reads.

supp_read_table <- tibble(
  Sample = colnames(countData),
  Assigned_reads = colSums(countData, na.rm = TRUE),
  Detected_features = colSums(countData > 0, na.rm = TRUE),
  Mean_count_per_detected_feature = apply(countData, 2, function(x) {
    mean(x[x > 0], na.rm = TRUE)
  }),
  Median_count_per_detected_feature = apply(countData, 2, function(x) {
    median(x[x > 0], na.rm = TRUE)
  })
)

# Add sample metadata
supp_read_table <- supp_read_table %>%
  left_join(
    colData %>%
      rownames_to_column("Sample"),
    by = "Sample"
  )

# -----------------------------
# 4. Optional: add featureCounts assignment summary
# -----------------------------
# This file is usually produced by featureCounts and often ends in ".summary".
# It contains Assigned, Unassigned_NoFeatures, Unassigned_Ambiguity, etc.
# Update the path below if your summary file has a different name.

featurecounts_summary_file <- 
  "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/Processed/A_Equina_Counts_unstranded.txt.summary"

if (file.exists(featurecounts_summary_file)) {
  
  fc_summary <- read.table(
    featurecounts_summary_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  
  fc_long <- fc_summary %>%
    pivot_longer(
      cols = -Status,
      names_to = "Sample_ID",
      values_to = "Read_count"
    ) %>%
    mutate(
      Sample_ID = basename(Sample_ID),
      Sample_ID = str_remove(Sample_ID, "\\.sorted.bam$"),
      Sample_ID = str_remove(Sample_ID, "\\.sorted.sam$")
    )
  
  fc_wide <- fc_long %>%
    pivot_wider(
      names_from = Status,
      values_from = Read_count
    ) %>%
    mutate(
      featureCounts_total_reads = rowSums(across(where(is.numeric)), na.rm = TRUE),
      featureCounts_assignment_rate_percent = 
        Assigned / featureCounts_total_reads * 100
    )
  
  fc_wide$Sample_ID <- as.integer(fc_wide$Sample_ID)
  
  supp_read_table <- supp_read_table %>%
    left_join(fc_wide, by = "Sample_ID")
}

# -----------------------------
# 5. Optional: add STAR alignment metrics
# -----------------------------
# This assumes each STAR output directory contains Log.final.out.
# Edit star_log_dir to the folder containing your STAR output folders or logs.

star_log_dir <- 
  "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/Processed/STAR_logs"

star_logs <- list.files(
  star_log_dir,
  pattern = "Log.final.out$",
  recursive = TRUE,
  full.names = TRUE
)

read_star_log <- function(log_file) {
  
  x <- readLines(log_file)
  
  get_metric <- function(pattern) {
    value <- x[str_detect(x, fixed(pattern))]
    if (length(value) == 0) return(NA_character_)
    str_split(value[1], "\\|", simplify = TRUE)[, 2] %>%
      str_trim()
  }
  
  tibble(
    Sample = basename(dirname(log_file)),
    Input_reads = as.numeric(get_metric("Number of input reads")),
    Uniquely_mapped_reads = as.numeric(get_metric("Uniquely mapped reads number")),
    Uniquely_mapped_percent = parse_number(get_metric("Uniquely mapped reads %")),
    Multi_mapped_percent = parse_number(get_metric("% of reads mapped to multiple loci")),
    Too_many_loci_percent = parse_number(get_metric("% of reads mapped to too many loci")),
    Unmapped_too_many_mismatches_percent = parse_number(get_metric("% of reads unmapped: too many mismatches")),
    Unmapped_too_short_percent = parse_number(get_metric("% of reads unmapped: too short")),
    Unmapped_other_percent = parse_number(get_metric("% of reads unmapped: other"))
  )
}

if (length(star_logs) > 0) {
  
  star_summary <- bind_rows(lapply(star_logs, read_star_log)) %>%
    mutate(
      Overall_alignment_rate_percent =
        Uniquely_mapped_percent + Multi_mapped_percent + Too_many_loci_percent
    )
  
  supp_read_table <- supp_read_table %>%
    left_join(star_summary, by = "Sample")
}

# -----------------------------
# 6. Clean and export table
# -----------------------------

supp_read_table <- supp_read_table %>%
  relocate(Sample, Diesel, Salinity, everything()) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write.csv(
  supp_read_table,
  "Supplementary_Table_Read_Summary.csv",
  row.names = FALSE
)

kable(
  supp_read_table,
  caption = "Supplementary table. Read count, assignment, and alignment summary for each sample."
)
