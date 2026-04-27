# Chapter 1: RNA-seq differential expression analysis

This directory contains the analysis code, data organisation, logs, figures, and results for Chapter 1 of the PhD thesis. Chapter 1 performs differential expression analysis on paired acute and paired primed RNA-seq experiments using DESeq2, followed by annotation, comparison of log2 fold-change patterns, and functional enrichment analysis.

The chapter directory is designed to keep raw data, processed data, scripts, results, figures, logs, and documentation separate so that the workflow is easier to follow and reproduce.

## Directory overview

```text
chapter_1/
├── README.md
├── data/
├── docs/
├── figures/
├── logs/
├── notebooks/
├── results/
└── scripts/
```

## `data/`

The `data/` directory contains input data and intermediate files used by the Chapter 1 RNA-seq workflow.

```text
data/
├── raw/
├── trimmed/
├── processed/
├── metadata/
├── external/
├── colData.tsv
├── featurecounts_alignment_summary.csv
├── hisat2_alignment_summary.csv
└── hisat2_merged_with_experimental_data_quotes.csv
```

### `data/raw/`

Expected to contain the original raw RNA-seq read files, typically FASTQ or FASTQ.GZ files. These files are large and are ignored by Git. The directory may contain a lightweight `README.md` placeholder in GitHub, but the sequencing files themselves should remain local or on HPC storage.

### `data/trimmed/`

Expected to contain trimmed RNA-seq reads generated after read-cleaning or adapter/quality trimming. These are also large sequencing files and are ignored by Git.

### `data/processed/`

Contains processed inputs and intermediate files used by downstream analysis. In the current workflow this directory is expected to contain BAM files used for read counting with featureCounts. It may also contain subdirectories such as:

```text
processed/
├── bam/
├── extracts/
└── index/
```

Expected contents include:

- aligned RNA-seq BAM files used as input to featureCounts
- featureCounts input or output intermediates
- local index files or alignment-derived files
- exon extract files or other processed annotation-derived files, where applicable

Large binary and alignment files such as `.bam`, `.bai`, and index files are ignored by Git.

### `data/processed/extracts/`

This directory is expected to contain processed exon extract files if they are generated locally. However, GTF and GFF3 genome annotation files are treated as external datasets and should be stored under `data/external/`.

### `data/external/`

Contains externally sourced reference and annotation files used in the Chapter 1 analysis. This includes genome annotation files such as `.gtf`, `.gff`, or `.gff3` files used by featureCounts.

The current featureCounts script expects annotation files to be located here:

```text
data/external/
```

and first looks for:

```text
data/external/A_Equina.gtf
```

If that file is not present, it will use the first `.gtf`, `.gff`, or `.gff3` file in `data/external/`.

### `data/metadata/`

Contains sample metadata and experiment-design files. This folder is intended for files describing sample IDs, treatment groups, pairing structure, experiment labels, and any other information needed to construct DESeq2 design formulae.

### Top-level files in `data/`

| File | Purpose |
|---|---|
| `colData.tsv` | Sample metadata used for DESeq2 analysis. This should contain sample-level information such as treatment, experiment, pairing/blocking variables, and other model covariates. |
| `featurecounts_alignment_summary.csv` | Summary of featureCounts assignment or counting statistics. |
| `hisat2_alignment_summary.csv` | Summary of HISAT2 alignment metrics. |
| `hisat2_merged_with_experimental_data_quotes.csv` | HISAT2 alignment summary merged with experimental metadata. |

## `scripts/`

The `scripts/` directory contains executable analysis scripts.

```text
scripts/
├── r/
└── unix/
```

### `scripts/r/`

Contains R and R Markdown scripts for the Chapter 1 RNA-seq analysis. Current scripts include:

| Script | Purpose |
|---|---|
| `Differential_Expression_Analysis_C1.Rmd` | Main DESeq2 differential expression analysis for the paired acute and primed experiments. |
| `Alignment_Efficiency_Analysis_C1.R` | Summarises and visualises alignment/counting efficiency metrics. |
| `Correlating_exp1_exp2_lfcs_C1.R` | Compares log2 fold-change estimates between the acute and primed experiments. |
| `Intersecting_Acute_DEA_results_with_and_without_outlier_C1.R` | Compares acute differential expression results with and without an outlier sample. |
| `Adding_real_gene_names_to_DEGs_lists_C1.R` | Adds gene or protein annotations to DEG result tables. |
| `Enrichment_Analyses_C1.Rmd` | Functional enrichment analysis for Chapter 1 gene sets. |
| `Enrichment_Analyses_new_blastp_C1.Rmd` | Updated enrichment workflow using BLASTP-derived annotation. |
| `Enrichment_Analyses_new_BLASTP_results.C1.Rmd` | Additional enrichment analysis using newer BLASTP results. |

### `scripts/unix/`

Contains shell and SLURM scripts used to run command-line analysis steps on the HPC.

Current script:

| Script | Purpose |
|---|---|
| `run_featurecounts_chapter_1.sh` | Runs featureCounts on Chapter 1 BAM files using annotation files from `data/external/`, writing the count matrix to `results/tables/featurecounts/`. |

The featureCounts script assumes it is submitted from the `chapter_1/` directory:

```bash
cd /uoa/scratch/users/r02hw22/repos/PhD---omics-analysis-portfolio/chapter_1
sbatch scripts/unix/run_featurecounts_chapter_1.sh
```

## `results/`

The `results/` directory contains generated analysis outputs from the Chapter 1 workflow.

```text
results/
├── tables/
├── enrichment/
├── gene_lists/
├── model_outputs/
└── figures/
```

### Top-level result files

The repository currently contains a mixture of final and intermediate result files at the top level of `results/`, including:

| File pattern | Purpose |
|---|---|
| `DEA_results_*.csv` | Differential expression result tables from DESeq2 analyses. |
| `DEresults_*.csv` | Additional DESeq2 result tables, including outlier-specific outputs. |
| `DEGs_1.txt`, `DEGs_2.txt` | Differentially expressed gene lists for the two Chapter 1 experiments. |
| `DE_results_with_*names*` | Differential expression results with added gene/protein annotation. |
| `GSEA_*` | Gene set enrichment analysis outputs. |
| `ORA_*` | Over-representation analysis outputs. |

### `results/tables/`

Intended for final tabular outputs used in the thesis or supplementary materials. The featureCounts script writes count matrices to:

```text
results/tables/featurecounts/
```

### `results/enrichment/`

Intended for functional enrichment outputs, including GO, GSEA, ORA, semantic clustering, and annotation-derived enrichment tables.

### `results/gene_lists/`

Intended for DEG lists, intersected gene lists, and other gene-set files used as input for downstream annotation or enrichment analysis.

### `results/model_outputs/`

Intended for model-derived outputs such as saved DESeq2 objects, transformed count matrices, model summaries, or other intermediate analysis objects. Large R objects such as `.RData` and `.rds` files are ignored by Git unless deliberately added.

### `results/figures/`

Contains generated figures that are outputs of the Chapter 1 analysis, including alignment/count summaries, PCA plots, volcano plots, heatmaps, Venn diagrams, enrichment plots, and log2 fold-change correlation plots.

Known subdirectories include:

```text
results/figures/
├── Heatmaps/
├── PCA plots/
├── Venn Diagrams/
└── Volcano Plots/
```

These folders contain analysis-generated plots. Some names still contain spaces and capital letters; these can be standardised later if required.

## `figures/`

The `figures/` directory is intended for thesis-organised figures rather than all generated analysis plots.

```text
figures/
├── exploratory/
├── main_text/
└── supplementary/
```

| Directory | Purpose |
|---|---|
| `figures/exploratory/` | Diagnostic and exploratory plots, such as PCA, QC, or preliminary visualisations. |
| `figures/main_text/` | Final figures intended for the main thesis chapter. |
| `figures/supplementary/` | Figures intended for supplementary materials. |

## `logs/`

The `logs/` directory stores scheduler and command-line job logs.

```text
logs/
├── outputs/
└── errors/
```

| Directory | Purpose |
|---|---|
| `logs/outputs/` | Standard output logs from SLURM jobs, such as featureCounts output logs. |
| `logs/errors/` | Standard error logs from SLURM jobs. |

Log files such as `.out` and `.err` are ignored by Git to avoid committing large or transient HPC job output.

## `docs/`

Contains written notes documenting Chapter 1 analysis decisions.

| File | Purpose |
|---|---|
| `methods_notes.md` | Notes on analysis methods and workflow decisions. |
| `contrast_notes.md` | Notes on DESeq2 contrasts and model comparisons. |
| `enrichment_notes.md` | Notes on functional enrichment analyses. |

## `notebooks/`

This directory is reserved for notebook-style analysis files or rendered exploratory notebooks. At present, most major Chapter 1 analysis notebooks are stored in `scripts/r/` as R Markdown files.

## Reproducibility notes

Large sequencing and alignment files are not expected to be tracked by GitHub. The `.gitignore` file is configured to ignore common large RNA-seq files and intermediate outputs, including:

- raw and trimmed FASTQ files
- BAM, BAI, SAM, and alignment files
- index files
- SLURM log files
- large R binary objects
- common RNA-seq and enrichment intermediate files

To reproduce the analysis, ensure the local HPC copy of `chapter_1/` contains the required ignored data files in the expected locations before running the scripts.

## Expected high-level workflow

1. Place raw RNA-seq reads in `data/raw/`.
2. Place trimmed reads in `data/trimmed/`.
3. Place aligned BAM files in `data/processed/` or `data/processed/bam/`.
4. Place GTF/GFF/GFF3 annotation files in `data/external/`.
5. Run featureCounts using `scripts/unix/run_featurecounts_chapter_1.sh`.
6. Use the resulting count matrix and `data/colData.tsv` for DESeq2 analysis in `scripts/r/Differential_Expression_Analysis_C1.Rmd`.
7. Run annotation, comparison, and enrichment scripts in `scripts/r/`.
8. Save generated tables under `results/` and thesis-ready plots under `figures/`.

## Directories needing confirmation

Some directory contents cannot be fully confirmed from GitHub because large files are intentionally ignored. Please verify locally that the following directories contain the expected files:

| Directory | Expected local contents |
|---|---|
| `data/raw/` | Raw FASTQ/FASTQ.GZ RNA-seq reads. |
| `data/trimmed/` | Trimmed FASTQ/FASTQ.GZ RNA-seq reads. |
| `data/processed/` | BAM files, featureCounts-relevant processed files, and possibly `extracts/` or `bam/` subdirectories. |
| `data/external/` | External annotation/reference files, especially GTF/GFF/GFF3 files. |
| `logs/outputs/` and `logs/errors/` | Local SLURM job logs, ignored by Git. |
