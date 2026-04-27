library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(patchwork)

# All lines on same plot - Use these figures ----------------------------------------------------
summary_df <- read.csv(
  "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/regionalised_aggregated_counts/combined_exp1_exp2_summary_df_for_plotting.csv"
)

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

# ---- colours
group_cols <- c(
  "control" = "#0047AB",  # strong royal blue
  "acute"   = "#E66100",  # vivid orange
  "naïve"   = "#009E73",  # bold green-teal
  "primed"  = "#CC007A"   # strong magenta
)

# ---- fixed y-range across all plots
y_max <- max(summary_df$max_count, na.rm = TRUE) + 1

# ---- make one stacked 3-context plot per region
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
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(aes(y = mean_count_log), linewidth = 1,
              alpha = 0.50) +
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
      breaks = seq(0, 100, 20)
    ) +
    scale_y_log10(limits = c(1, y_max)) +
    labs(
      title = region_labels[[region_name]],
      x = "% methylation",
      y = expression(log[10]("Count"))
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

# ---- combine into 4 x 2 layout with shared legend 
combined_plot <- wrap_plots(region_plots, ncol = 4, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot

ggsave(
  filename = file.path("C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/regionalised_aggregated_counts/exp1_and_exp2_stacked_methylation_percents.png"),
  plot = combined_plot,
  width = 14,
  height = 10,
  dpi = 600
)

# subsetting experiments from the new plotting data ----
summary_df <- read.csv(
  "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/regionalised_aggregated_counts/combined_exp1_exp2_summary_df_for_plotting.csv"
)

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

# ---- colours
group_cols <- c(
  "control" = "#0047AB",  # strong royal blue
  "acute"   = "#d91c3f",  # vivid orange
  "naïve"   = "#ffb300",  # bold green-teal
  "primed"  = "#009E73"   # strong magenta
)

# ---- fixed y-range across all plots
y_max <- max(summary_df$max_count, na.rm = TRUE) + 1

## exp1 ----
region_plots <- lapply(regions, function(region_name) {
  region_df <- summary_df %>%
    filter(experiment == "exp1", region == region_name) %>%
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
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(aes(y = mean_count_log), linewidth = 1,
              alpha = 0.50) +
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
      breaks = seq(0, 100, 20)
    ) +
    scale_y_log10(limits = c(1, y_max)) +
    labs(
      title = region_labels[[region_name]],
      x = "% methylation",
      y = expression(log[10]("Count"))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.spacing.y = grid::unit(0.2, "lines"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
})

# ---- combine into 4 x 2 layout with shared legend 
combined_plot_exp1 <- wrap_plots(region_plots, ncol = 4, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot_exp1

ggsave(
  filename = file.path("C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/regionalised_aggregated_counts/exp1_stacked_methylation_percents.png"),
  plot = combined_plot_exp1,
  width = 14,
  height = 10,
  dpi = 600
)

## exp2 ----
region_plots <- lapply(regions, function(region_name) {
  region_df <- summary_df %>%
    filter(experiment == "exp2", region == region_name) %>%
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
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(aes(y = mean_count_log), linewidth = 1,
              alpha = 0.50) +
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
      breaks = seq(0, 100, 20)
    ) +
    scale_y_log10(limits = c(1, y_max)) +
    labs(
      title = region_labels[[region_name]],
      x = "% methylation",
      y = expression(log[10]("Count"))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.spacing.y = grid::unit(0.2, "lines"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
})

# ---- combine into 4 x 2 layout with shared legend 
combined_plot_exp2 <- wrap_plots(region_plots, ncol = 4, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot_exp2

ggsave(
  filename = file.path("C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/regionalised_aggregated_counts/exp2_stacked_methylation_percents.png"),
  plot = combined_plot_exp2,
  width = 14,
  height = 10,
  dpi = 600
)

## Combined ----
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
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(aes(y = mean_count_log), linewidth = 1,
              alpha = 0.50) +
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
      breaks = seq(0, 100, 20)
    ) +
    scale_y_log10(limits = c(1, y_max)) +
    labs(
      title = region_labels[[region_name]],
      x = "% methylation",
      y = expression(log[10]("Count"))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.spacing.y = grid::unit(0.2, "lines"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
})

# ---- combine into 4 x 2 layout with shared legend 
combined_plot <- wrap_plots(region_plots, ncol = 4, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot

ggsave(
  filename = file.path("C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/regionalised_aggregated_counts/exp1_and_exp2_stacked_methylation_percents.png"),
  plot = combined_plot,
  width = 14,
  height = 10,
  dpi = 600
)


# Code to generate data from aggregated methylation counts files: ----
library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(patchwork)

## ---- paths and settings ----
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

## ---- sample metadata from your experiment setup ----
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

## ---- helper: robust file reader ----
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

## ---- helper: find likely column names ----
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

## ---- helper: extract % methylation ----
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

## ---- helper: extract file prefix from filename ----
extract_file_prefix <- function(file_path, region) {
  basename(file_path) %>%
    str_remove(paste0("_", region, "_aggregated_methylation_counts\\.txt$"))
}

## ---- read and combine all experiments + regions + contexts ----
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

## ---- add sample metadata and clean up ----
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

## ---- binning ----
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

## ---- count per sample within each region/context/group ----
counts_df <- binned_df %>%
  count(region, context, experiment, group, sample_id, bin_mid, name = "count") %>%
  group_by(region, context, experiment, group, sample_id) %>%
  complete(bin_mid = bin_mids, fill = list(count = 0)) %>%
  ungroup()

## ---- summarise counts across samples within each group ----
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

## ---- save plotting dataframes ----
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

