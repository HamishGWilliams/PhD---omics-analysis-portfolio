!/usr/bin/env Rscript
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

# ---- Setup ----
base_path <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results"
new_base_path <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures"

experiments <- c("exp1", "exp2")
contexts <- c("CpG", "CHG", "CHH")

setwd(base_path)
genome <- read_gff3("../../Data2/combined_annotations.gff3")
regions <- unique(genome$type)

dir.create(new_base_path, recursive = TRUE, showWarnings = FALSE)

copy_log <- data.frame(
  experiment = character(),
  context = character(),
  region = character(),
  sample_id = character(),
  source_file = character(),
  dest_file = character(),
  copied = logical(),
  stringsAsFactors = FALSE
)

for (experiment in experiments) {
  for (context in contexts) {
    for (region in regions) {

      exp_context_region_path <- file.path(base_path, experiment, context, region)

      if (experiment == "exp1") {
        file.list <- c(
          paste0("3_",  region, "_aggregated_methylation_counts.txt"),
          paste0("9_",  region, "_aggregated_methylation_counts.txt"),
          paste0("13_", region, "_aggregated_methylation_counts.txt"),
          paste0("15_", region, "_aggregated_methylation_counts.txt"),
          paste0("16_", region, "_aggregated_methylation_counts.txt"),
          paste0("17_", region, "_aggregated_methylation_counts.txt"),
          paste0("22_", region, "_aggregated_methylation_counts.txt"),
          paste0("27_", region, "_aggregated_methylation_counts.txt")
        )
        sample.id <- c("3", "9", "13", "15", "16", "17", "22", "27")

      } else if (experiment == "exp2") {
        file.list <- c(
          paste0("6_",     region, "_aggregated_methylation_counts.txt"),
          paste0("7_",     region, "_aggregated_methylation_counts.txt"),
          paste0("8_",     region, "_aggregated_methylation_counts.txt"),
          paste0("10_",    region, "_aggregated_methylation_counts.txt"),
          paste0("14_",    region, "_aggregated_methylation_counts.txt"),
          paste0("20_",    region, "_aggregated_methylation_counts.txt"),
          paste0("21_",    region, "_aggregated_methylation_counts.txt"),
          paste0("29_11_", region, "_aggregated_methylation_counts.txt"),
          paste0("32_18_", region, "_aggregated_methylation_counts.txt"),
          paste0("36_2_",  region, "_aggregated_methylation_counts.txt")
        )
        sample.id <- c("6", "7", "8", "10", "14", "20", "21", "11", "18", "2")
      }

      dest_dir <- file.path(new_base_path, experiment, context, region)
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

      for (i in seq_along(file.list)) {
        source_file <- file.path(exp_context_region_path, file.list[i])
        dest_file <- file.path(dest_dir, file.list[i])

        copied <- FALSE
        if (file.exists(source_file)) {
          copied <- file.copy(from = source_file, to = dest_file, overwrite = TRUE)
        } else {
          warning("Missing file: ", source_file)
        }

        copy_log <- rbind(
          copy_log,
          data.frame(
            experiment = experiment,
            context = context,
            region = region,
            sample_id = sample.id[i],
            source_file = source_file,
            dest_file = dest_file,
            copied = copied,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

write.csv(copy_log, file.path(new_base_path, "copied_files_manifest.csv"), row.names = FALSE)

cat("Done. Files copied to:\n", new_base_path, "\n")


# Making the first histogram

library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

# ---- paths ----
data_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures/exp1/CpG/gene"

files <- list.files(
  path = data_dir,
  pattern = "_gene_aggregated_methylation_counts\\.txt$",
  full.names = TRUE
)

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  # first try tab-delimited
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  # if everything came in as one column, try generic whitespace-delimited
  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  # if still one column, try base read.table
  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  return(df)
}

# ---- helper to find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  return(NA_character_)
}

# ---- helper to extract % methylation ----
extract_percent_methylation <- function(df) {

  # direct percentage columns
  direct_col <- find_first_matching_col(df, c(
    "percent",
    "percentage",
    "pct",
    "prop",
    "proportion",
    "weighted.*meth",
    "meth.*percent",
    "methylation"
  ))

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))

    # keep only if at least some values are numeric
    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  # methylated count columns
  meth_col <- find_first_matching_col(df, c(
    "^m$",
    "^mc$",
    "^meth$",
    "methylated",
    "num.*meth",
    "count.*meth",
    "c_count",
    "numc"
  ))

  # unmethylated count columns
  unmeth_col <- find_first_matching_col(df, c(
    "^u$",
    "^uc$",
    "^unmeth$",
    "unmethylated",
    "num.*unmeth",
    "count.*unmeth",
    "t_count",
    "numt"
  ))

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  # methylated + total coverage
  total_col <- find_first_matching_col(df, c(
    "^n$",
    "^cov$",
    "coverage",
    "total",
    "depth",
    "num.*total"
  ))

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- inspect first file before combining ----
test_df <- read_meth_file(files[1])
cat("Detected columns in first file:\n")
print(colnames(test_df))
print(head(test_df))

# ---- combine all files ----
plot_df <- map_dfr(files, function(f) {
  df <- read_meth_file(f)

  tibble(
    sample_id = str_extract(basename(f), "^\\d+"),
    percent_methylation = extract_percent_methylation(df)
  )
})

plot_df <- plot_df %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100)

# ---- plot ----
p <- ggplot(plot_df, aes(x = percent_methylation)) +
  geom_histogram(binwidth = 5, boundary = 0, closed = "left") +
  facet_wrap(~ sample_id, ncol = 3) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "exp1 | CpG | gene",
    x = "% methylation per region",
    y = "Number of gene regions"
  ) +
  theme_bw(base_size = 12)

# save plot to file
ggsave(
  filename = file.path(data_dir, "exp1_CpG_gene_methylation_histograms.png"),
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)

# single plot instead of a facet
plot_df$sample_id <- factor(plot_df$sample_id)

p_overlay <- ggplot(plot_df, aes(x = percent_methylation, fill = sample_id)) +
  geom_histogram(
    binwidth = 5,
    boundary = 0,
    closed = "left",
    position = "identity",
    alpha = 0.4
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "exp1 | CpG | gene",
    x = "% methylation per region",
    y = "Number of gene regions",
    fill = "Sample ID"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(data_dir, "exp1_CpG_gene_methylation_histogram_overlay.png"),
  plot = p_overlay,
  width = 10,
  height = 7,
  dpi = 300
)

# change to a line graph ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

plot_df$sample_id <- factor(plot_df$sample_id)

binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2   # 2.5, 7.5, 12.5, ...

# assign each row to a bin index, then map that index to the midpoint
binned_df <- plot_df %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# count regions per sample per bin
counts_df <- binned_df %>%
  count(sample_id, bin_mid, name = "count")

# keep zero-count bins for each sample
counts_df <- counts_df %>%
  complete(
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# summarise across samples
summary_df <- counts_df %>%
  group_by(bin_mid) %>%
  summarise(
    mean_count = mean(count),
    min_count = min(count),
    max_count = max(count),
    .groups = "drop"
  )

# plot
p_summary <- ggplot(summary_df, aes(x = bin_mid, y = mean_count)) +
  geom_ribbon(aes(ymin = min_count, ymax = max_count), alpha = 0.25) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "exp1 | CpG | gene",
    x = "% methylation per region",
    y = "Mean number of gene regions across samples"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(data_dir, "exp1_CpG_gene_mean_histogram_range.png"),
  plot = p_summary,
  width = 10,
  height = 7,
  dpi = 300
)

# change y-axis to proportion -------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)

plot_df$sample_id <- factor(plot_df$sample_id)

binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

# assign each region to a bin
binned_df <- plot_df %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# count regions per sample per bin
counts_df <- binned_df %>%
  count(sample_id, bin_mid, name = "count") %>%
  complete(
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# convert counts to proportions within each sample
props_df <- counts_df %>%
  group_by(sample_id) %>%
  mutate(
    total_regions = sum(count),
    proportion = count / total_regions
  ) %>%
  ungroup()

# summarise proportions across samples
summary_df <- props_df %>%
  group_by(bin_mid) %>%
  summarise(
    mean_prop = mean(proportion),
    min_prop  = min(proportion),
    max_prop  = max(proportion),
    .groups = "drop"
  )

# plot
p_summary_prop <- ggplot(summary_df, aes(x = bin_mid, y = mean_prop)) +
  geom_ribbon(aes(ymin = min_prop, ymax = max_prop), alpha = 0.25) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "exp1 | CpG | gene",
    x = "% methylation per region",
    y = "Proportion of gene regions"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(data_dir, "exp1_CpG_gene_mean_proportion_histogram_range.png"),
  plot = p_summary_prop,
  width = 10,
  height = 7,
  dpi = 300
)




# Adding all Cytosine contexts "CpG", "CHG", "CHH" ------------------------------------------------

library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# ---- paths and settings ----
base_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures/exp1"
contexts <- c("CpG", "CHG", "CHH")
region <- "gene"

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  df
}

# ---- helper: find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  NA_character_
}

# ---- helper: extract % methylation ----
extract_percent_methylation <- function(df) {

  direct_col <- find_first_matching_col(df, c(
    "percent",
    "percentage",
    "pct",
    "prop",
    "proportion",
    "weighted.*meth",
    "meth.*percent",
    "methylation"
  ))

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))

    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  meth_col <- find_first_matching_col(df, c(
    "^m$",
    "^mc$",
    "^meth$",
    "methylated",
    "num.*meth",
    "count.*meth",
    "c_count",
    "numc"
  ))

  unmeth_col <- find_first_matching_col(df, c(
    "^u$",
    "^uc$",
    "^unmeth$",
    "unmethylated",
    "num.*unmeth",
    "count.*unmeth",
    "t_count",
    "numt"
  ))

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  total_col <- find_first_matching_col(df, c(
    "^n$",
    "^cov$",
    "coverage",
    "total",
    "depth",
    "num.*total"
  ))

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- read and combine all contexts ----
plot_df_all <- map_dfr(contexts, function(context) {

  data_dir <- file.path(base_dir, context, region)

  files <- list.files(
    path = data_dir,
    pattern = "_gene_aggregated_methylation_counts\\.txt$",
    full.names = TRUE
  )

  map_dfr(files, function(f) {
    df <- read_meth_file(f)

    tibble(
      context = context,
      sample_id = str_extract(basename(f), "^\\d+"),
      percent_methylation = extract_percent_methylation(df)
    )
  })
})

# clean up
plot_df_all <- plot_df_all %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100) %>%
  mutate(
    context = factor(context, levels = c("CpG", "CHG", "CHH")),
    sample_id = factor(sample_id)
  )

# ---- binning ----
binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

binned_df <- plot_df_all %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# ---- count per sample within each context ----
counts_df <- binned_df %>%
  count(context, sample_id, bin_mid, name = "count") %>%
  complete(
    context,
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# ---- convert counts to within-sample proportions ----
props_df <- counts_df %>%
  group_by(context, sample_id) %>%
  mutate(
    total_regions = sum(count),
    proportion = count / total_regions
  ) %>%
  ungroup()

# ---- summarise across samples within each context ----
summary_df <- props_df %>%
  group_by(context, bin_mid) %>%
  summarise(
    mean_prop = mean(proportion),
    min_prop  = min(proportion),
    max_prop  = max(proportion),
    .groups = "drop"
  )

# ---- stacked context plot ----
p_summary_prop <- ggplot(summary_df, aes(x = bin_mid, y = mean_prop)) +
  geom_ribbon(aes(ymin = min_prop, ymax = max_prop), alpha = 0.25) +
  geom_line(linewidth = 1) +
  facet_grid(context ~ .) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "exp1 | gene | methylation distribution by context",
    x = "% methylation per region",
    y = "Proportion of gene regions"
  ) +
  theme_bw(base_size = 12)

# ---- save ----
save_dir <- file.path(base_dir, "stacked_context_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(save_dir, "exp1_gene_contexts_stacked_mean_proportion_histogram_range.png"),
  plot = p_summary_prop,
  width = 10,
  height = 12,
  dpi = 300
)

# Switch y-axis back to true values instead of proportions --------------------------------------


library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# ---- paths and settings ----
base_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures/exp1"
contexts <- c("CpG", "CHG", "CHH")
region <- "gene"

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  df
}

# ---- helper: find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  NA_character_
}

# ---- helper: extract % methylation ----
extract_percent_methylation <- function(df) {

  direct_col <- find_first_matching_col(df, c(
    "percent",
    "percentage",
    "pct",
    "prop",
    "proportion",
    "weighted.*meth",
    "meth.*percent",
    "methylation"
  ))

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))

    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  meth_col <- find_first_matching_col(df, c(
    "^m$",
    "^mc$",
    "^meth$",
    "methylated",
    "num.*meth",
    "count.*meth",
    "c_count",
    "numc"
  ))

  unmeth_col <- find_first_matching_col(df, c(
    "^u$",
    "^uc$",
    "^unmeth$",
    "unmethylated",
    "num.*unmeth",
    "count.*unmeth",
    "t_count",
    "numt"
  ))

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  total_col <- find_first_matching_col(df, c(
    "^n$",
    "^cov$",
    "coverage",
    "total",
    "depth",
    "num.*total"
  ))

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- read and combine all contexts ----
plot_df_all <- map_dfr(contexts, function(context) {

  data_dir <- file.path(base_dir, context, region)

  files <- list.files(
    path = data_dir,
    pattern = "_gene_aggregated_methylation_counts\\.txt$",
    full.names = TRUE
  )

  map_dfr(files, function(f) {
    df <- read_meth_file(f)

    tibble(
      context = context,
      sample_id = str_extract(basename(f), "^\\d+"),
      percent_methylation = extract_percent_methylation(df)
    )
  })
})

# clean up
plot_df_all <- plot_df_all %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100) %>%
  mutate(
    context = factor(context, levels = c("CpG", "CHG", "CHH")),
    sample_id = factor(sample_id)
  )

# ---- binning ----
binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

binned_df <- plot_df_all %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# ---- count per sample within each context ----
counts_df <- binned_df %>%
  count(context, sample_id, bin_mid, name = "count") %>%
  complete(
    context,
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# ---- summarise raw counts across samples ----
summary_df <- counts_df %>%
  group_by(context, bin_mid) %>%
  summarise(
    mean_count = mean(count),
    min_count  = min(count),
    max_count  = max(count),
    .groups = "drop"
  )

# ---- stacked context plot ----
p_summary_counts <- ggplot(
  summary_df,
  aes(x = bin_mid, y = mean_count, colour = context, fill = context)
) +
  geom_ribbon(aes(ymin = min_count, ymax = max_count), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1) +
  facet_grid(context ~ .) +
  scale_color_manual(values = c(
    "CpG" = "#0072B2",
    "CHG" = "#D55E00",
    "CHH" = "#009E73"
  )) +
  scale_fill_manual(values = c(
    "CpG" = "#0072B2",
    "CHG" = "#D55E00",
    "CHH" = "#009E73"
  )) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "exp1 | gene | methylation distribution by context",
    x = "% methylation per region",
    y = "Number of gene regions"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
# ---- save ----
save_dir <- file.path(base_dir, "stacked_context_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(save_dir, "exp1_gene_contexts_stacked_mean_count_histogram_range.png"),
  plot = p_summary_counts,
  width = 10,
  height = 12,
  dpi = 300
)


# All region types on one plot ------------------------------------------------------------------


library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# ---- paths and settings ----
base_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures/exp1"
contexts <- c("CpG", "CHG", "CHH")

# specify only the region types you want
regions <- c(
  "dispersed_repeat",
  "downstream_region",
  "exon",
  "five_prime_UTR",
  "gene",
  "ncRNA_gene",
  "promoter",
  "three_prime_UTR"
)

# optional safety check: keep only regions that actually exist
regions <- regions[regions %in% list.dirs(
  file.path(base_dir, "CpG"),
  recursive = FALSE,
  full.names = FALSE
)]

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  df
}

# ---- helper: find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  NA_character_
}

# ---- helper: extract % methylation ----
extract_percent_methylation <- function(df) {

  direct_col <- find_first_matching_col(df, c(
    "percent",
    "percentage",
    "pct",
    "prop",
    "proportion",
    "weighted.*meth",
    "meth.*percent",
    "methylation"
  ))

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))

    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  meth_col <- find_first_matching_col(df, c(
    "^m$",
    "^mc$",
    "^meth$",
    "methylated",
    "num.*meth",
    "count.*meth",
    "c_count",
    "numc"
  ))

  unmeth_col <- find_first_matching_col(df, c(
    "^u$",
    "^uc$",
    "^unmeth$",
    "unmethylated",
    "num.*unmeth",
    "count.*unmeth",
    "t_count",
    "numt"
  ))

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  total_col <- find_first_matching_col(df, c(
    "^n$",
    "^cov$",
    "coverage",
    "total",
    "depth",
    "num.*total"
  ))

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- read and combine all regions + contexts ----
plot_df_all <- map_dfr(regions, function(region) {
  map_dfr(contexts, function(context) {

    data_dir <- file.path(base_dir, context, region)

    files <- list.files(
      path = data_dir,
      pattern = paste0("_", region, "_aggregated_methylation_counts\\.txt$"),
      full.names = TRUE
    )

    map_dfr(files, function(f) {
      df <- read_meth_file(f)

      tibble(
        region = region,
        context = context,
        sample_id = str_extract(basename(f), "^\\d+"),
        percent_methylation = extract_percent_methylation(df)
      )
    })
  })
})

# ---- clean up ----
plot_df_all <- plot_df_all %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100) %>%
  mutate(
    region = factor(region, levels = regions),
    context = factor(context, levels = c("CpG", "CHG", "CHH")),
    sample_id = factor(sample_id)
  )

# ---- binning ----
binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

binned_df <- plot_df_all %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# ---- count per sample within each region/context ----
counts_df <- binned_df %>%
  count(region, context, sample_id, bin_mid, name = "count") %>%
  complete(
    region,
    context,
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# ---- summarise raw counts across samples ----
summary_df <- counts_df %>%
  group_by(region, context, bin_mid) %>%
  summarise(
    mean_count = mean(count),
    min_count  = min(count),
    max_count  = max(count),
    .groups = "drop"
  )

# ---- plot: all region types in 4 x 2 layout ----
p_summary_counts <- ggplot(
  summary_df,
  aes(x = bin_mid, y = mean_count, colour = context, fill = context)
) +
  geom_ribbon(aes(ymin = min_count, ymax = max_count), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ region, ncol = 4, nrow = 2, scales = "free_y") +
  scale_color_manual(values = c(
    "CpG" = "#0072B2",
    "CHG" = "#D55E00",
    "CHH" = "#009E73"
  )) +
  scale_fill_manual(values = c(
    "CpG" = "#0072B2",
    "CHG" = "#D55E00",
    "CHH" = "#009E73"
  )) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "exp1 | methylation distribution across region types",
    x = "% methylation per region",
    y = "Number of regions",
    colour = "Context",
    fill = "Context"
  ) +
  theme_bw(base_size = 12)

# ---- save ----
save_dir <- file.path(base_dir, "all_region_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(save_dir, "exp1_all_regions_4x2_mean_count_histogram_range.png"),
  plot = p_summary_counts,
  width = 16,
  height = 10,
  dpi = 300
)

# Plots with triple stacked panels ---------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(patchwork)

# ---- paths and settings ----
base_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures/exp1"

contexts <- c("CpG", "CHG", "CHH")

regions <- c(
  "dispersed_repeat",
  "downstream_region",
  "exon",
  "five_prime_UTR",
  "gene",
  "ncRNA_gene",
  "promoter",
  "three_prime_UTR"
)

region_labels <- c(
  "dispersed_repeat"   = "Dispersed repeat",
  "downstream_region"  = "Downstream region",
  "exon"               = "Exon",
  "five_prime_UTR"     = "5' UTR",
  "gene"               = "Gene",
  "ncRNA_gene"         = "ncRNA gene",
  "promoter"           = "Promoter",
  "three_prime_UTR"    = "3' UTR"
)

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  df
}

# ---- helper: find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  NA_character_
}

# ---- helper: extract % methylation ----
extract_percent_methylation <- function(df) {
  direct_col <- find_first_matching_col(
    df,
    c(
      "percent",
      "percentage",
      "pct",
      "prop",
      "proportion",
      "weighted.*meth",
      "meth.*percent",
      "methylation"
    )
  )

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))

    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  meth_col <- find_first_matching_col(
    df,
    c(
      "^m$",
      "^mc$",
      "^meth$",
      "methylated",
      "num.*meth",
      "count.*meth",
      "c_count",
      "numc"
    )
  )

  unmeth_col <- find_first_matching_col(
    df,
    c(
      "^u$",
      "^uc$",
      "^unmeth$",
      "unmethylated",
      "num.*unmeth",
      "count.*unmeth",
      "t_count",
      "numt"
    )
  )

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  total_col <- find_first_matching_col(
    df,
    c("^n$", "^cov$", "coverage", "total", "depth", "num.*total")
  )

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- read and combine all regions + contexts ----
plot_df_all <- map_dfr(regions, function(region) {
  map_dfr(contexts, function(context) {
    data_dir <- file.path(base_dir, context, region)

    files <- list.files(
      path = data_dir,
      pattern = paste0("_", region, "_aggregated_methylation_counts\\.txt$"),
      full.names = TRUE
    )

    map_dfr(files, function(f) {
      df <- read_meth_file(f)

      tibble(
        region = region,
        context = context,
        sample_id = str_extract(basename(f), "^\\d+"),
        percent_methylation = extract_percent_methylation(df)
      )
    })
  })
})

# ---- clean up ----
plot_df_all <- plot_df_all %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100) %>%
  mutate(
    region = factor(region, levels = regions),
    context = factor(context, levels = c("CpG", "CHG", "CHH")),
    sample_id = factor(sample_id)
  )

# ---- binning ----
binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

binned_df <- plot_df_all %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# ---- count per sample within each region/context ----
counts_df <- binned_df %>%
  count(region, context, sample_id, bin_mid, name = "count") %>%
  complete(
    region,
    context,
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# ---- summarise raw counts across samples ----
summary_df <- counts_df %>%
  group_by(region, context, bin_mid) %>%
  summarise(
    mean_count = mean(count),
    min_count = min(count),
    max_count = max(count),
    .groups = "drop"
  )

# ---- save plotting dataframes ----
save_dir <- file.path(base_dir, "all_region_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# final dataframe used for plotting
write.csv(
  summary_df,
  file.path(save_dir, "exp1_summary_df_for_plotting.csv"),
  row.names = FALSE
)

# optional: save the underlying per-sample binned counts too
write.csv(
  counts_df,
  file.path(save_dir, "exp1_counts_df_for_plotting.csv"),
  row.names = FALSE
)

# optional: save the raw combined methylation values
write.csv(
  plot_df_all,
  file.path(save_dir, "exp1_plot_df_all_raw.csv"),
  row.names = FALSE
)

# ---- colours ----
context_cols <- c(
  "CpG" = "#0072B2",
  "CHG" = "#D55E00",
  "CHH" = "#009E73"
)

# ---- make one stacked 3-context plot per region ----
region_plots <- lapply(regions, function(region_name) {
  region_df <- summary_df %>%
    filter(region == region_name)

  ggplot(
    region_df,
    aes(x = bin_mid, y = mean_count, colour = context, fill = context)
  ) +
    geom_ribbon(
      aes(ymin = min_count, ymax = max_count),
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    facet_grid(
      rows = vars(context),
      scales = "free_y",
      switch = "y"
    ) +
    scale_color_manual(values = context_cols) +
    scale_fill_manual(values = context_cols) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    labs(
      title = region_labels[[region_name]],
      x = "% methylation per region",
      y = "Number of regions"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      panel.spacing.y = grid::unit(0.2, "lines")
    )
})

# ---- combine into 4 x 2 layout ----
combined_plot <- wrap_plots(region_plots, ncol = 4, nrow = 2)

# ---- save ----
save_dir <- file.path(base_dir, "all_region_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(save_dir, "exp1_all_regions_nested_4x2_stacked_contexts.png"),
  plot = combined_plot,
  width = 20,
  height = 14,
  dpi = 300
)



# exp2 plot now ------------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(patchwork)

# ---- paths and settings ----
base_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures/exp2"

contexts <- c("CpG", "CHG", "CHH")

regions <- c(
  "dispersed_repeat",
  "downstream_region",
  "exon",
  "five_prime_UTR",
  "gene",
  "ncRNA_gene",
  "promoter",
  "three_prime_UTR"
)

region_labels <- c(
  "dispersed_repeat"   = "Dispersed repeat",
  "downstream_region"  = "Downstream region",
  "exon"               = "Exon",
  "five_prime_UTR"     = "5' UTR",
  "gene"               = "Gene",
  "ncRNA_gene"         = "ncRNA gene",
  "promoter"           = "Promoter",
  "three_prime_UTR"    = "3' UTR"
)

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  df
}

# ---- helper: find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  NA_character_
}

# ---- helper: extract % methylation ----
extract_percent_methylation <- function(df) {
  direct_col <- find_first_matching_col(
    df,
    c(
      "percent",
      "percentage",
      "pct",
      "prop",
      "proportion",
      "weighted.*meth",
      "meth.*percent",
      "methylation"
    )
  )

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))

    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  meth_col <- find_first_matching_col(
    df,
    c(
      "^m$",
      "^mc$",
      "^meth$",
      "methylated",
      "num.*meth",
      "count.*meth",
      "c_count",
      "numc"
    )
  )

  unmeth_col <- find_first_matching_col(
    df,
    c(
      "^u$",
      "^uc$",
      "^unmeth$",
      "unmethylated",
      "num.*unmeth",
      "count.*unmeth",
      "t_count",
      "numt"
    )
  )

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  total_col <- find_first_matching_col(
    df,
    c("^n$", "^cov$", "coverage", "total", "depth", "num.*total")
  )

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- read and combine all regions + contexts ----
plot_df_all <- map_dfr(regions, function(region) {
  map_dfr(contexts, function(context) {
    data_dir <- file.path(base_dir, context, region)

    files <- list.files(
      path = data_dir,
      pattern = paste0("_", region, "_aggregated_methylation_counts\\.txt$"),
      full.names = TRUE
    )

    map_dfr(files, function(f) {
      df <- read_meth_file(f)

      tibble(
        region = region,
        context = context,
        sample_id = str_extract(basename(f), "^\\d+"),
        percent_methylation = extract_percent_methylation(df)
      )
    })
  })
})

# ---- clean up ----
plot_df_all <- plot_df_all %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100) %>%
  mutate(
    region = factor(region, levels = regions),
    context = factor(context, levels = c("CpG", "CHG", "CHH")),
    sample_id = factor(sample_id)
  )

# ---- binning ----
binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

binned_df <- plot_df_all %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# ---- count per sample within each region/context ----
counts_df <- binned_df %>%
  count(region, context, sample_id, bin_mid, name = "count") %>%
  complete(
    region,
    context,
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# ---- summarise raw counts across samples ----
summary_df <- counts_df %>%
  group_by(region, context, bin_mid) %>%
  summarise(
    mean_count = mean(count),
    min_count = min(count),
    max_count = max(count),
    .groups = "drop"
  )

# ---- save plotting dataframes ----
save_dir <- file.path(base_dir, "all_region_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# final dataframe used for plotting
write.csv(
  summary_df,
  file.path(save_dir, "exp2_summary_df_for_plotting.csv"),
  row.names = FALSE
)

# optional: save the underlying per-sample binned counts too
write.csv(
  counts_df,
  file.path(save_dir, "exp2_counts_df_for_plotting.csv"),
  row.names = FALSE
)

# optional: save the raw combined methylation values
write.csv(
  plot_df_all,
  file.path(save_dir, "exp2_plot_df_all_raw.csv"),
  row.names = FALSE
)

# ---- colours ----
context_cols <- c(
  "CpG" = "#0072B2",
  "CHG" = "#D55E00",
  "CHH" = "#009E73"
)

# ---- make one stacked 3-context plot per region ----
region_plots <- lapply(regions, function(region_name) {
  region_df <- summary_df %>%
    filter(region == region_name)

  ggplot(
    region_df,
    aes(x = bin_mid, y = mean_count, colour = context, fill = context)
  ) +
    geom_ribbon(
      aes(ymin = min_count, ymax = max_count),
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    facet_grid(
      rows = vars(context),
      scales = "free_y",
      switch = "y"
    ) +
    scale_color_manual(values = context_cols) +
    scale_fill_manual(values = context_cols) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    labs(
      title = region_labels[[region_name]],
      x = "% methylation per region",
      y = "Number of regions"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      panel.spacing.y = grid::unit(0.2, "lines")
    )
})

# ---- combine into 4 x 2 layout ----
combined_plot <- wrap_plots(region_plots, ncol = 4, nrow = 2)

# ---- save ----
save_dir <- file.path(base_dir, "all_region_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(save_dir, "exp1_all_regions_nested_4x2_stacked_contexts.png"),
  plot = combined_plot,
  width = 20,
  height = 14,
  dpi = 300
)



# Separate lines for control and treatment conditions ----------------------------------------------


library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(patchwork)

# ---- paths and settings ----
base_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures"
experiment <- "exp2"

contexts <- c("CpG", "CHG", "CHH")

regions <- c(
  "dispersed_repeat",
  "downstream_region",
  "exon",
  "five_prime_UTR",
  "gene",
  "ncRNA_gene",
  "promoter",
  "three_prime_UTR"
)

region_labels <- c(
  "dispersed_repeat"  = "Dispersed repeat",
  "downstream_region" = "Downstream region",
  "exon"              = "Exon",
  "five_prime_UTR"    = "5' UTR",
  "gene"              = "Gene",
  "ncRNA_gene"        = "ncRNA gene",
  "promoter"          = "Promoter",
  "three_prime_UTR"   = "3' UTR"
)

# ---- sample metadata from your original experiment setup ----
if (experiment == "exp1") {
  sample_info <- tibble(
    file_prefix = c("9", "13", "15", "16", "22", "27"),
    sample_id   = c("9", "13", "15", "16", "22", "27"),
    treatment   = c(0, 0, 1, 0, 1, 1)
  )
} else if (experiment == "exp2") {
  sample_info <- tibble(
    file_prefix = c("6", "7", "10", "14", "20", "21", "29_11", "36_2"),
    sample_id   = c("6", "7", "10", "14", "20", "21", "11", "2"),
    treatment   = c(0, 0, 0, 1, 1, 1, 0, 1)
  )
} else {
  stop("experiment must be either 'exp1' or 'exp2'")
}

sample_info <- sample_info %>%
  mutate(
    group = factor(
      if_else(treatment == 0, "Control", "Treatment"),
      levels = c("Control", "Treatment")
    )
  )

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  df
}

# ---- helper: find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  NA_character_
}

# ---- helper: extract % methylation ----
extract_percent_methylation <- function(df) {
  direct_col <- find_first_matching_col(
    df,
    c(
      "percent",
      "percentage",
      "pct",
      "prop",
      "proportion",
      "weighted.*meth",
      "meth.*percent",
      "methylation"
    )
  )

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))

    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  meth_col <- find_first_matching_col(
    df,
    c(
      "^m$",
      "^mc$",
      "^meth$",
      "methylated",
      "num.*meth",
      "count.*meth",
      "c_count",
      "numc"
    )
  )

  unmeth_col <- find_first_matching_col(
    df,
    c(
      "^u$",
      "^uc$",
      "^unmeth$",
      "unmethylated",
      "num.*unmeth",
      "count.*unmeth",
      "t_count",
      "numt"
    )
  )

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  total_col <- find_first_matching_col(
    df,
    c("^n$", "^cov$", "coverage", "total", "depth", "num.*total")
  )

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- helper: extract file prefix from filename ----
extract_file_prefix <- function(file_path, region) {
  basename(file_path) %>%
    str_remove(paste0("_", region, "_aggregated_methylation_counts\\.txt$"))
}

# ---- read and combine all regions + contexts ----
plot_df_all <- map_dfr(regions, function(region) {
  map_dfr(contexts, function(context) {
    data_dir <- file.path(base_dir, experiment, context, region)

    files <- list.files(
      path = data_dir,
      pattern = paste0("_", region, "_aggregated_methylation_counts\\.txt$"),
      full.names = TRUE
    )

    map_dfr(files, function(f) {
      file_prefix <- extract_file_prefix(f, region)

      if (!file_prefix %in% sample_info$file_prefix) {
        return(tibble())
      }

      df <- read_meth_file(f)

      tibble(
        region = region,
        context = context,
        file_prefix = file_prefix,
        percent_methylation = extract_percent_methylation(df)
      )
    })
  })
})

# ---- add sample metadata and clean up ----
plot_df_all <- plot_df_all %>%
  left_join(sample_info, by = "file_prefix") %>%
  filter(!is.na(sample_id)) %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100) %>%
  mutate(
    region = factor(region, levels = regions),
    context = factor(context, levels = c("CpG", "CHG", "CHH")),
    sample_id = factor(sample_id, levels = unique(sample_info$sample_id)),
    group = factor(group, levels = c("Control", "Treatment"))
  )

# ---- binning ----
binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

binned_df <- plot_df_all %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# ---- count per sample within each region/context/group ----
counts_df <- binned_df %>%
  count(region, context, group, sample_id, bin_mid, name = "count") %>%
  complete(
    region,
    context,
    group,
    sample_id,
    bin_mid = bin_mids,
    fill = list(count = 0)
  )

# ---- summarise counts across samples within each group ----
summary_df <- counts_df %>%
  group_by(region, context, group, bin_mid) %>%
  summarise(
    mean_count = mean(count),
    min_count = min(count),
    max_count = max(count),
    .groups = "drop"
  )

# ---- save plotting dataframes ----
save_dir <- file.path(base_dir, experiment, "all_region_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  summary_df,
  file.path(save_dir, paste0(experiment, "_summary_df_for_plotting.csv")),
  row.names = FALSE
)

write.csv(
  counts_df,
  file.path(save_dir, paste0(experiment, "_counts_df_for_plotting.csv")),
  row.names = FALSE
)

write.csv(
  plot_df_all,
  file.path(save_dir, paste0(experiment, "_plot_df_all_raw.csv")),
  row.names = FALSE
)

saveRDS(
  summary_df,
  file.path(save_dir, paste0(experiment, "_summary_df_for_plotting.rds"))
)

saveRDS(
  counts_df,
  file.path(save_dir, paste0(experiment, "_counts_df_for_plotting.rds"))
)

saveRDS(
  plot_df_all,
  file.path(save_dir, paste0(experiment, "_plot_df_all_raw.rds"))
)

# ---- colours ----
group_cols <- c(
  "Control" = "#0072B2",
  "Treatment" = "#D55E00"
)

# ---- fixed y-range across all plots ----
y_max <- max(summary_df$max_count, na.rm = TRUE) + 1

# ---- make one stacked 3-context plot per region ----
region_plots <- lapply(regions, function(region_name) {
  region_df <- summary_df %>%
    filter(region == region_name) %>%
    mutate(
      mean_count_log = mean_count + 1,
      min_count_log  = min_count + 1,
      max_count_log  = max_count + 1
    )

  ggplot(
    region_df,
    aes(
      x = bin_mid,
      y = mean_count_log,
      colour = group,
      fill = group,
      group = group
    )
  ) +
    geom_ribbon(
      aes(ymin = min_count_log, ymax = max_count_log),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    facet_grid(
      rows = vars(context),
      scales = "fixed",
      switch = "y"
    ) +
    scale_color_manual(values = group_cols, name = "Group") +
    scale_fill_manual(values = group_cols, name = "Group") +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    scale_y_log10(
      limits = c(1, y_max)
    ) +
    labs(
      title = region_labels[[region_name]],
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      strip.text.y.left = element_blank(),
      panel.spacing.y = grid::unit(0.2, "lines"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
})

# ---- combine into 4 x 2 layout with shared legend ----
combined_plot <- wrap_plots(region_plots, ncol = 4, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ---- shared axis labels ----
shared_x <- grid::textGrob(
  "% methylation per region",
  gp = grid::gpar(fontsize = 12)
)

shared_y <- grid::textGrob(
  expression(log[10]("Number of regions" + 1)),
  rot = 90,
  gp = grid::gpar(fontsize = 12)
)

combined_plot_with_axes <-
  (wrap_elements(shared_y) + combined_plot) /
  wrap_elements(shared_x) +
  plot_layout(
    widths = c(0.04, 1),
    heights = c(1, 0.05)
  )

# ---- save plots ----
ggsave(
  filename = file.path(
    save_dir,
    paste0(experiment, "_all_regions_nested_4x2_stacked_contexts_control_vs_treatment.png")
  ),
  plot = combined_plot_with_axes,
  width = 20,
  height = 14,
  dpi = 300
)


# All lines on same plot ----------------------------------------------------------------------


library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(patchwork)

# ---- paths and settings ----
base_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results_for_figures"

experiments <- c("exp1", "exp2")
contexts <- c("CpG", "CHG", "CHH")

regions <- c(
  "dispersed_repeat",
  "downstream_region",
  "exon",
  "five_prime_UTR",
  "gene",
  "ncRNA_gene",
  "promoter",
  "three_prime_UTR"
)

region_labels <- c(
  "dispersed_repeat" = "Dispersed repeat",
  "downstream_region" = "Downstream region",
  "exon" = "Exon",
  "five_prime_UTR" = "5' UTR",
  "gene" = "Gene",
  "ncRNA_gene" = "ncRNA gene",
  "promoter" = "Promoter",
  "three_prime_UTR" = "3' UTR"
)

# ---- sample metadata from your experiment setup ----
sample_info <- bind_rows(
  tibble(
    experiment = "exp1",
    file_prefix = c("9", "13", "15", "16", "22", "27"),
    sample_id = c("9", "13", "15", "16", "22", "27"),
    treatment = c(0, 0, 1, 0, 1, 1)
  ),
  tibble(
    experiment = "exp2",
    file_prefix = c("6", "7", "10", "14", "20", "21", "29_11", "36_2"),
    sample_id = c("6", "7", "10", "14", "20", "21", "11", "2"),
    treatment = c(0, 0, 0, 1, 1, 1, 0, 1)
  )
) %>%
  mutate(
    group = case_when(
      experiment == "exp1" & treatment == 0 ~ "control",
      experiment == "exp1" & treatment == 1 ~ "acute",
      experiment == "exp2" & treatment == 0 ~ "naïve",
      experiment == "exp2" & treatment == 1 ~ "primed",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("control", "acute", "naïve", "primed"))
  )

# ---- helper: robust file reader ----
read_meth_file <- function(f) {
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE))

  if (ncol(df) == 1) {
    df <- suppressMessages(read_table(f, show_col_types = FALSE))
  }

  if (ncol(df) == 1) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df <- as_tibble(df)
  }

  df
}

# ---- helper: find likely column names ----
find_first_matching_col <- function(df, patterns) {
  nms <- names(df)
  low <- tolower(nms)

  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }

  NA_character_
}

# ---- helper: extract % methylation ----
extract_percent_methylation <- function(df) {
  direct_col <- find_first_matching_col(
    df,
    c(
      "percent", "percentage", "pct", "prop", "proportion",
      "weighted.*meth", "meth.*percent", "methylation"
    )
  )

  if (!is.na(direct_col)) {
    vals <- suppressWarnings(as.numeric(df[[direct_col]]))
    if (sum(!is.na(vals)) > 0) {
      if (all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
        vals <- vals * 100
      }
      return(vals)
    }
  }

  meth_col <- find_first_matching_col(
    df,
    c(
      "^m$", "^mc$", "^meth$", "methylated", "num.*meth",
      "count.*meth", "c_count", "numc"
    )
  )

  unmeth_col <- find_first_matching_col(
    df,
    c(
      "^u$", "^uc$", "^unmeth$", "unmethylated", "num.*unmeth",
      "count.*unmeth", "t_count", "numt"
    )
  )

  if (!is.na(meth_col) && !is.na(unmeth_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    unmeth <- suppressWarnings(as.numeric(df[[unmeth_col]]))
    return(100 * meth / (meth + unmeth))
  }

  total_col <- find_first_matching_col(
    df,
    c("^n$", "^cov$", "coverage", "total", "depth", "num.*total")
  )

  if (!is.na(meth_col) && !is.na(total_col)) {
    meth <- suppressWarnings(as.numeric(df[[meth_col]]))
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    return(100 * meth / total)
  }

  stop(
    paste0(
      "Could not identify columns needed to calculate % methylation.\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  )
}

# ---- helper: extract file prefix from filename ----
extract_file_prefix <- function(file_path, region) {
  basename(file_path) %>%
    str_remove(paste0("_", region, "_aggregated_methylation_counts\\.txt$"))
}

# ---- read and combine all experiments + regions + contexts ----
plot_df_all <- map_dfr(experiments, function(experiment_name) {
  map_dfr(regions, function(region) {
    map_dfr(contexts, function(context) {
      data_dir <- file.path(base_dir, experiment_name, context, region)

      files <- list.files(
        path = data_dir,
        pattern = paste0("_", region, "_aggregated_methylation_counts\\.txt$"),
        full.names = TRUE
      )

      map_dfr(files, function(f) {
        file_prefix <- extract_file_prefix(f, region)

        if (!file_prefix %in% sample_info$file_prefix[sample_info$experiment == experiment_name]) {
          return(tibble())
        }

        df <- read_meth_file(f)

        tibble(
          experiment = experiment_name,
          region = region,
          context = context,
          file_prefix = file_prefix,
          percent_methylation = extract_percent_methylation(df)
        )
      })
    })
  })
})

# ---- add sample metadata and clean up ----
plot_df_all <- plot_df_all %>%
  left_join(sample_info, by = c("experiment", "file_prefix")) %>%
  filter(!is.na(sample_id), !is.na(group)) %>%
  filter(is.finite(percent_methylation)) %>%
  filter(percent_methylation >= 0, percent_methylation <= 100) %>%
  mutate(
    region = factor(region, levels = regions),
    context = factor(context, levels = c("CpG", "CHG", "CHH")),
    group = factor(group, levels = c("control", "acute", "naïve", "primed"))
  )

# ---- binning ----
binwidth <- 5
breaks <- seq(0, 100, by = binwidth)
bin_mids <- breaks[-1] - binwidth / 2

binned_df <- plot_df_all %>%
  mutate(
    bin_index = cut(
      percent_methylation,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    ),
    bin_mid = bin_mids[bin_index]
  )

# ---- count per sample within each region/context/group ----
counts_df <- binned_df %>%
  count(region, context, experiment, group, sample_id, bin_mid, name = "count") %>%
  group_by(region, context, experiment, group, sample_id) %>%
  complete(bin_mid = bin_mids, fill = list(count = 0)) %>%
  ungroup()

# ---- summarise counts across samples within each group ----
summary_df <- counts_df %>%
  group_by(experiment, region, context, group, bin_mid) %>%
  summarise(
    mean_count = mean(count),
    min_count = min(count),
    max_count = max(count),
    .groups = "drop"
  ) %>%
  mutate(
    experiment = factor(experiment, levels = c("exp1", "exp2")),
    group = factor(group, levels = c("control", "acute", "naïve", "primed"))
  )

# ---- save plotting dataframes ----
save_dir <- file.path(base_dir, "combined_exp1_exp2_plots")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  summary_df,
  file.path(save_dir, "combined_exp1_exp2_summary_df_for_plotting.csv"),
  row.names = FALSE
)
write.csv(
  counts_df,
  file.path(save_dir, "combined_exp1_exp2_counts_df_for_plotting.csv"),
  row.names = FALSE
)
write.csv(
  plot_df_all,
  file.path(save_dir, "combined_exp1_exp2_plot_df_all_raw.csv"),
  row.names = FALSE
)

saveRDS(
  summary_df,
  file.path(save_dir, "combined_exp1_exp2_summary_df_for_plotting.rds")
)
saveRDS(
  counts_df,
  file.path(save_dir, "combined_exp1_exp2_counts_df_for_plotting.rds")
)
saveRDS(
  plot_df_all,
  file.path(save_dir, "combined_exp1_exp2_plot_df_all_raw.rds")
)

# ---- colours ----
group_cols <- c(
  "control" = "#0072B2",
  "acute" = "#D55E00",
  "naïve" = "#009E73",
  "primed" = "#CC79A7"
)

# ---- fixed y-range across all plots ----
y_max <- max(summary_df$max_count, na.rm = TRUE) + 1

# ---- make one stacked 3-context plot per region ----
region_plots <- lapply(regions, function(region_name) {
  region_df <- summary_df %>%
    filter(region == region_name) %>%
    mutate(
      mean_count_log = mean_count + 1,
      min_count_log = min_count + 1,
      max_count_log = max_count + 1,
      ribbon_ymin = ifelse(min_count_log <= 1, 1.02, min_count_log),
      ribbon_ymax = max_count_log
    )

  ggplot(region_df, aes(
    x = bin_mid,
    colour = group,
    fill = group,
    group = group
  )) +
    geom_ribbon(
      aes(ymin = ribbon_ymin, ymax = ribbon_ymax),
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(aes(y = mean_count_log), linewidth = 1) +
    facet_grid(
      rows = vars(context),
      scales = "fixed",
      switch = "y"
    ) +
    scale_color_manual(
      values = group_cols,
      breaks = c("control", "acute", "naïve", "primed"),
      labels = c("control", "acute", "naïve", "primed"),
      name = "Group"
    ) +
    scale_fill_manual(
      values = group_cols,
      breaks = c("control", "acute", "naïve", "primed"),
      labels = c("control", "acute", "naïve", "primed"),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    scale_y_log10(limits = c(1, y_max)) +
    labs(
      title = region_labels[[region_name]],
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      strip.text.y.left = element_blank(),
      panel.spacing.y = grid::unit(0.2, "lines"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
})

# ---- combine into 4 x 2 layout with shared legend ----
combined_plot <- wrap_plots(region_plots, ncol = 4, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ---- shared axis labels ----
shared_x <- grid::textGrob(
  "% methylation per region",
  gp = grid::gpar(fontsize = 12)
)

shared_y <- grid::textGrob(
  expression(log[10]("Number of regions" + 1)),
  rot = 90,
  gp = grid::gpar(fontsize = 12)
)

combined_plot_with_axes <- (wrap_elements(shared_y) + combined_plot) /
  wrap_elements(shared_x) +
  plot_layout(
    widths = c(0.04, 1),
    heights = c(1, 0.05)
  )

# ---- save plots ----
ggsave(
  filename = file.path(
    save_dir,
    "combined_exp1_exp2_all_regions_nested_4x2_stacked_contexts.png"
  ),
  plot = combined_plot_with_axes,
  width = 20,
  height = 14,
  dpi = 300
)

ggsave(
  filename = file.path(
    save_dir,
    "combined_exp1_exp2_all_regions_nested_4x2_stacked_contexts.pdf"
  ),
  plot = combined_plot_with_axes,
  width = 20,
  height = 14
)






