# PhD Omics Analysis Portfolio

This repository contains the analysis code and supporting material for the three data chapters of my PhD thesis. The project brings together transcriptomic and DNA methylation analyses used to investigate molecular responses across paired acute, paired priming, and multi-stressor exposure experiments.

The repository is intended to provide a transparent and reproducible record of the analytical workflows used throughout the thesis, including quality control, differential expression analysis, differential methylation analysis, functional enrichment, data integration, and figure generation.

## Project overview

The thesis is structured around three data chapters:

| Chapter | Focus | Main analysis type |
|---|---|---|
| Chapter 1 | Gene expression responses in paired acute and paired primed experiments | Differential expression analysis using DESeq2 |
| Chapter 2 | DNA methylation responses in the same paired acute and paired primed experiments | Differential methylation analysis and comparison with gene expression |
| Chapter 3 | Gene expression responses to combined diesel and salinity exposure | Differential expression analysis using DESeq2 |

## Repository structure

```text
.
├── Chapter 1/
│   └── Scripts/
│       └── R/
│           ├── Differential_Expression_Analysis_C1.Rmd
│           ├── Enrichment_Analyses_C1.Rmd
│           └── Additional scripts for QC, DEG processing, and visualisation
│
├── Chapter 2/
│   └── Scripts/
│       └── R/
│           ├── Methylation_Analysis_Final_C2.Rmd
│           ├── DMR functional enrichment analyses.Rmd
│           └── Additional scripts for methylation summaries and visualisation
│
├── Chapter 3/
│   └── Code/
│       ├── Differential_Expression_Analysis_C3.Rmd
│       ├── Enrichment_analyses_C3.Rmd
│       └── Additional scripts for PCA, upset plots, and model comparisons
│
└── README.md
```

## Chapter 1: Differential expression analysis of paired acute and primed experiments

Chapter 1 investigates transcriptomic responses using RNA-seq data from two paired experimental designs:

1. A paired acute exposure experiment
2. A paired primed exposure experiment

Differential expression analysis is performed using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), with models designed to account for the paired structure of the experiments. The analyses include:

- Import and preparation of gene count matrices
- Sample metadata formatting
- Quality control and exploratory visualisation
- Principal component analysis
- Differential expression testing
- Identification of differentially expressed genes
- Comparison of log2 fold-change patterns between experimental contrasts
- Functional enrichment analysis
- Generation of thesis figures and supporting outputs

The main analysis file for this chapter is:

```text
Chapter 1/Scripts/R/Differential_Expression_Analysis_C1.Rmd
```

## Chapter 2: Differential methylation analysis and integration with expression

Chapter 2 investigates DNA methylation responses in the same paired acute and primed experimental contexts analysed in Chapter 1.

The analysis focuses on identifying differentially methylated regions or loci, annotating methylation changes, and comparing methylation patterns with gene expression responses. This chapter therefore provides an integrated view of transcriptional and epigenetic responses across the paired experimental designs.

The analyses include:

- Processing and summarising methylation data
- Differential methylation analysis
- Annotation of methylated regions or loci
- Functional enrichment of methylation-associated genes
- Comparison of methylation and expression patterns
- Visualisation of methylation percentages and genomic feature associations
- Generation of thesis figures and supporting tables

The main analysis file for this chapter is:

```text
Chapter 2/Scripts/R/Methylation_Analysis_Final_C2.Rmd
```

## Chapter 3: Differential expression analysis of multi-stressor exposure

Chapter 3 analyses a separate multi-stressor experiment investigating transcriptomic responses to diesel and salinity exposure.

This chapter uses RNA-seq differential expression analysis to examine how gene expression changes under single and combined stressor conditions. The analysis is designed to identify genes and pathways associated with diesel exposure, salinity exposure, and their combined effects.

The analyses include:

- Import and preparation of RNA-seq count data
- Sample-level quality control
- Principal component analysis
- Differential expression analysis using DESeq2
- Model comparison and contrast-specific testing
- Identification of treatment-responsive genes
- Functional enrichment analysis
- Visualisation of shared and treatment-specific responses
- Generation of thesis figures and supporting outputs

The main analysis file for this chapter is:

```text
Chapter 3/Code/Differential_Expression_Analysis_C3.Rmd
```

## Software and dependencies

Most analyses are written in R and R Markdown. Key packages used across the repository include:

- `DESeq2` for RNA-seq differential expression analysis
- `tidyverse` for data wrangling and visualisation
- `ggplot2` for figure generation
- `pheatmap` or related packages for heatmaps
- `clusterProfiler`, `topGO`, or equivalent tools for enrichment analysis
- Additional chapter-specific packages for methylation analysis, annotation, and plotting

Package requirements may vary between chapters. See the individual R or R Markdown scripts for chapter-specific libraries.

## Reproducibility notes

The repository is organised by thesis chapter, with each chapter containing the scripts required to reproduce the corresponding analyses and figures.

To reproduce an analysis:

1. Clone the repository:

```bash
git clone https://github.com/HamishGWilliams/PhD---omics-analysis-portfolio.git
```

2. Open the relevant chapter directory.

3. Review the corresponding R or R Markdown script.

4. Check that required input files are available in the expected locations.

5. Run the scripts in the order indicated by filenames, comments, or the relevant chapter workflow.

Large raw sequencing files are not necessarily stored in this repository. Where applicable, scripts assume that processed count matrices, methylation files, metadata tables, or annotation files are available locally.

## Outputs

The analyses generate outputs including:

- Differential expression result tables
- Differential methylation result tables
- Gene lists
- Annotated gene or region tables
- Enrichment results
- PCA plots
- Volcano plots
- Upset plots
- Methylation summary plots
- Thesis-ready figures and supplementary outputs

## Purpose of the repository

This repository serves three main purposes:

1. To archive the computational analyses used in my PhD thesis
2. To provide a transparent record of the decisions made during analysis
3. To act as a portfolio of applied omics workflows across transcriptomic and epigenomic datasets

## Citation

If using or adapting code from this repository, please cite the associated thesis or contact the repository author.

```text
Williams, H. G. PhD thesis. [Thesis title and institution to be added]
```

## Author

**Hamish G. Williams**

This repository forms part of my PhD research portfolio and contains analysis code associated with transcriptomic and DNA methylation analyses across multiple experimental chapters.
