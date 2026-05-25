#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Combined volcano plot from existing methylKit DMA outputs
# ------------------------------------------------------------
# This script does NOT rerun methylKit.
# It reads existing differential methylation result tables
# and creates one faceted volcano plot containing all
# experiments and cytosine contexts.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

project_dir <- "/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2"

results_dir <- file.path(project_dir, "results")

model_output_dir <- file.path(
  results_dir,
  "model_outputs",
  "base_level_dma_covariates_no_outlier_pairs"
)

volcano_figure_dir <- file.path(
  results_dir,
  "figures",
  "volcano_base_level_covariates_no_outlier_pairs"
)

dir.create(volcano_figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------

experiments <- c("exp1", "exp2")
contexts <- c("CpG", "CHG", "CHH")

# ------------------------------------------------------------
# Volcano plot settings
# ------------------------------------------------------------

p_adj_threshold <- 0.1

x_limits <- c(-40, 40)
y_limits <- c(0, 6)

x_breaks <- seq(-40, 40, by = 10)
y_breaks <- seq(0, 10, by = 1)

point_alpha <- 0.6
point_size <- 1

output_width <- 24
output_height <- 30
output_units <- "cm"
output_dpi <- 300

# ------------------------------------------------------------
# Read one DMA table
# ------------------------------------------------------------

read_dma_table <- function(experiment, context) {

  input_file <- file.path(
    model_output_dir,
    paste0(
      "Actinia_DMBs_d5_",
      experiment,
      "_",
      context,
      "_covariates_no_outlier_pairs_raw.txt"
    )
  )

  # Fallback to annotated file if raw file is missing
  if (!file.exists(input_file)) {
    fallback_file <- file.path(
      model_output_dir,
      paste0(
        "act_diff_df_ann_",
        experiment,
        "_",
        context,
        "_covariates_no_outlier_pairs.txt"
      )
    )

    if (file.exists(fallback_file)) {
      warning(
        "Raw file not found. Using annotated file instead: ",
        fallback_file
      )
      input_file <- fallback_file
    } else {
      warning(
        "No input file found for ",
        experiment,
        " / ",
        context,
        ". Skipping."
      )
      return(NULL)
    }
  }

  message("Reading: ", input_file)

  dma_data <- read_tsv(
    input_file,
    col_types = cols(),
    show_col_types = FALSE
  )

  required_cols <- c("meth.diff", "pvalue", "p_fdr")
  missing_cols <- setdiff(required_cols, names(dma_data))

  if (length(missing_cols) > 0) {
    stop(
      "Input file is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nFile: ",
      input_file
    )
  }

  dma_data %>%
    mutate(
      experiment = experiment,
      context = context
    )
}

# ------------------------------------------------------------
# Combine all DMA tables
# ------------------------------------------------------------

all_volcano_data <- list()

for (experiment in experiments) {
  for (context in contexts) {
    table_i <- tryCatch(
      {
        read_dma_table(
          experiment = experiment,
          context = context
        )
      },
      error = function(e) {
        message(
          "Error reading table for ",
          experiment,
          " / ",
          context,
          ": ",
          e$message
        )
        NULL
      }
    )

    if (!is.null(table_i)) {
      all_volcano_data[[paste(experiment, context, sep = "_")]] <- table_i
    }
  }
}

if (length(all_volcano_data) == 0) {
  stop("No DMA tables could be read. No volcano plot generated.")
}

volcano_data <- bind_rows(all_volcano_data)

# ------------------------------------------------------------
# Prepare volcano plotting data
# ------------------------------------------------------------

volcano_data <- volcano_data %>%
  mutate(
    meth.diff = as.numeric(meth.diff),
    pvalue = as.numeric(pvalue),
    p_fdr = as.numeric(p_fdr),

    pvalue_plot = pmax(pvalue, .Machine$double.xmin),
    neg_log10_pvalue = -log10(pvalue_plot),

    padj_status = case_when(
      p_fdr <= p_adj_threshold ~ "padj <= 0.1",
      TRUE ~ "padj > 0.1"
    ),

    padj_status = factor(
      padj_status,
      levels = c("padj > 0.1", "padj <= 0.1")
    ),

    experiment = factor(
      experiment,
      levels = c("exp1", "exp2"),
      labels = c("Acute", "Primed")
    ),

    context = factor(
      context,
      levels = c("CpG", "CHG", "CHH")
    ),

    facet_label = paste(context, experiment, sep = " - ")
  ) %>%
  filter(
    is.finite(meth.diff),
    is.finite(neg_log10_pvalue),
    is.finite(p_fdr)
  )

if (nrow(volcano_data) == 0) {
  stop("No finite rows available for volcano plot.")
}

volcano_data$facet_label <- factor(
  volcano_data$facet_label,
  levels = c(
    "CpG - Acute",
    "CpG - Primed",
    "CHG - Acute",
    "CHG - Primed",
    "CHH - Acute",
    "CHH - Primed"
  )
)
# ------------------------------------------------------------
# Combined faceted volcano plot
# ------------------------------------------------------------

combined_volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = meth.diff,
    y = neg_log10_pvalue,
    colour = padj_status
  )
) +
  geom_point(
    alpha = point_alpha,
    size = point_size
  ) +
  geom_vline(
    xintercept = 0,
    colour = "black",
    linetype = "dotted"
  ) +
  facet_wrap(
    ~ facet_label,
    ncol = 2
  ) +
  scale_colour_manual(
    values = c(
      "padj > 0.1" = "grey70",
      "padj <= 0.1" = "red"
    ),
    labels = c(
      "padj > 0.1",
      "padj <= 0.1"
    ),
    drop = FALSE
  ) +
  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits
  ) +
  scale_x_continuous(
    breaks = x_breaks
  ) +
  scale_y_continuous(
    breaks = y_breaks
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(
      face = "bold",
      margin = margin(0, 20, 0, 0),
      size = rel(1.1),
      colour = "black"
    ),
    axis.title.x = element_text(
      hjust = 0.5,
      face = "bold",
      margin = margin(20, 0, 0, 0),
      size = rel(1.1),
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = rel(1.1)
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  ) +
  labs(
    x = "Differential methylation %",
    y = expression("-log"[10] * " raw p-value"),
    colour = "Adjusted p-value"
  )
# ------------------------------------------------------------
# Save plot
# ------------------------------------------------------------

output_file <- file.path(
  volcano_figure_dir,
  "combined_volcano_plot_all_experiments_contexts_covariates_no_outlier_pairs.png"
)

ggsave(
  filename = output_file,
  plot = combined_volcano_plot,
  width = output_width,
  height = output_height,
  units = output_units,
  dpi = output_dpi
)

message("Saved combined volcano plot: ", output_file)
message("Combined volcano plot complete.")