#!/bin/bash
#SBATCH --job-name=plot_mbias
#SBATCH --mem=16G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/logs/errors/%x_%j.err

module load r/4.2.2

EXTRACT_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/data/processed/methylation_extraction"
PLOT_DIR="/uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_2/results/figures/m_bias"

mkdir -p "$PLOT_DIR"

shopt -s nullglob

for sample_dir in "$EXTRACT_DIR"/*/ ; do
    sample=$(basename "$sample_dir")

    echo "Processing M-bias plots for sample: $sample"
    echo "Input directory: $sample_dir"

    mbias_files=("$sample_dir"/*.M-bias.txt)

    if [ "${#mbias_files[@]}" -eq 0 ]; then
        echo "No M-bias file found for sample: $sample; skipping."
        echo
        continue
    fi

    for mbias_file in "${mbias_files[@]}" ; do
        mbias_base=$(basename "$mbias_file" .M-bias.txt)
        sample_plot_dir="$PLOT_DIR/$sample"
        tmp_dir="$sample_dir/m_bias_split_tmp"

        mkdir -p "$sample_plot_dir"
        rm -rf "$tmp_dir"
        mkdir -p "$tmp_dir"

        echo "Splitting M-bias file: $mbias_file"

        awk -v out="$tmp_dir/file" '
            BEGIN { RS="================"; ORS="" }
            {
                f = out NR
                print $0 > f
                close(f)
            }
        ' "$mbias_file"

        echo "Generating plot for: $sample"

        Rscript --vanilla - "$tmp_dir" "$sample" "$mbias_base" "$sample_plot_dir" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

tmp_dir <- args[1]
sample <- args[2]
mbias_base <- args[3]
sample_plot_dir <- args[4]

suppressPackageStartupMessages(library(ggplot2))

read_mbias_section <- function(file, context, read) {
  if (!file.exists(file)) {
    warning("Missing expected M-bias section file: ", file)
    return(NULL)
  }

  dat <- tryCatch(
    read.table(
      file,
      header = FALSE,
      skip = 2,
      fill = TRUE,
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      warning("Could not read section file: ", file)
      return(NULL)
    }
  )

  if (is.null(dat) || nrow(dat) == 0 || ncol(dat) < 5) {
    warning("Section file was empty or malformed: ", file)
    return(NULL)
  }

  dat <- dat[, 1:5]
  colnames(dat) <- c("position", "meth", "un", "pc_meth", "coverage")

  dat$position <- suppressWarnings(as.numeric(dat$position))
  dat$meth <- suppressWarnings(as.numeric(dat$meth))
  dat$un <- suppressWarnings(as.numeric(dat$un))
  dat$pc_meth <- suppressWarnings(as.numeric(dat$pc_meth))
  dat$coverage <- suppressWarnings(as.numeric(dat$coverage))

  dat <- dat[!is.na(dat$position), ]

  if (nrow(dat) == 0) {
    warning("No numeric position rows found in section file: ", file)
    return(NULL)
  }

  dat$context <- context
  dat$read <- read

  dat
}

sections <- list(
  read_mbias_section(file.path(tmp_dir, "file2"), "CpG", "Read 1"),
  read_mbias_section(file.path(tmp_dir, "file3"), "CHG", "Read 1"),
  read_mbias_section(file.path(tmp_dir, "file4"), "CHH", "Read 1"),
  read_mbias_section(file.path(tmp_dir, "file5"), "CpG", "Read 2"),
  read_mbias_section(file.path(tmp_dir, "file6"), "CHG", "Read 2"),
  read_mbias_section(file.path(tmp_dir, "file7"), "CHH", "Read 2")
)

all_data <- do.call(rbind, sections[!vapply(sections, is.null, logical(1))])

if (is.null(all_data) || nrow(all_data) == 0) {
  stop("No usable M-bias data found for sample: ", sample)
}

all_data$context <- factor(all_data$context, levels = c("CpG", "CHG", "CHH"))
all_data$read <- factor(all_data$read, levels = c("Read 1", "Read 2"))

read_meth_pc <- ggplot(
  all_data,
  aes(
    x = position,
    y = pc_meth,
    group = context,
    colour = context
  )
) +
  geom_line(linewidth = 0.4) +
  theme_classic() +
  labs(
    title = paste("M-bias methylation profile:", sample),
    x = "Position in read",
    y = "Percentage methylated",
    colour = "Context"
  ) +
  facet_wrap(~read, ncol = 1)

out_png <- file.path(
  sample_plot_dir,
  paste0(mbias_base, "_m_bias_percent_methylated_by_position.png")
)

png(
  filename = out_png,
  width = 8,
  height = 6,
  units = "in",
  res = 300
)

print(read_meth_pc)
dev.off()

message("Saved plot: ", out_png)
RSCRIPT

        rm -rf "$tmp_dir"

        echo "Finished plot for: $sample"
        echo
    done
done

echo "All M-bias plots complete."