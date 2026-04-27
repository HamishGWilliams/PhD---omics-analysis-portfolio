# ---- load packages ----

# This function will install bioconductor packages also without needing to do the whole BiocManager::install(...), which is nice
check_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      print(paste(pkg, "is not installed, installing package..."))
      
      # Check if the package requires BiocManager
      if (pkg %in% c("gage", "RUVSeq", "gprofiler2","DESeq2","PCAtools", "ComplexHeatmap")) {
        if (!require("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
          BiocManager::install(version = "3.17")
        }
        BiocManager::install(pkg, force = TRUE)
      } else {
        install.packages(pkg)
      }
      
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
      print(paste(pkg, "is installed and loaded"))
    }
  }
}

# Generate a package list: 
pkg_list<-c("tidyverse", #manipulation
            "DESeq2", # main differential expression analysis tool
            "tibble", # tables
            "stats", # stats
            "ggplot2", # plotting
            "pheatmap", # plotting
            "dplyr", # manipulation 
            "tidyr", # manipulation
            "ggfortify", #plotting
            "corrplot", # plotting
            "gprofiler2", #
            "knitr", # manipulation
            "ggforce", # pltting
            "gridExtra", # plotting
            "ggrepel", # plotting
            "VennDiagram", # plotting
            "ggVennDiagram", # plotting
            "PCAtools", # plotting
            "cowplot", # plotting
            "apeglm", # LFC shrinking method option
            "EnhancedVolcano", # plotting
            "ashr", # LFC shrinking method option
            "ComplexHeatmap",
            "RColorBrewer"
) # plotting

# Execute function:
check_packages(pkg_list)


# ---- Set-up ----

# We will need results with both non-shrunk and shrunk lfc values
# therefore we need to repeat the DEA analysis and save objects with the raw and shrunk lfc values

# Import Count Data
# .rmd file is in .../Chapter 1/Scripts/
counts <- as.matrix(read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Data/A_Equina_Counts_redo.tsv", header = T, sep = '\t'))
counts_exp1 <- subset(counts[,c(1:5,6:12,27:33)])
# Remove the 'width' column
countData_exp1_with_outlier <- as.matrix(subset(counts_exp1, select = c(-Chr,-Start,-End,-Strand,-Length)))
# Remove outliers in Count Data
countData_exp1_wo_outlier <- as.matrix(subset(countData_exp1_with_outlier[,c(1:5, 7:9, 11:14)]))

# Here we tell R to make sure that it reads our data as numeric, as I had issues with this previously
{apply(countData_exp1_with_outlier, 2, as.numeric)
  sapply(countData_exp1_with_outlier, as.numeric)
  class(countData_exp1_with_outlier) <- "numeric"
  storage.mode(countData_exp1_with_outlier) <- "numeric"}

{apply(countData_exp1_wo_outlier, 2, as.numeric)
  sapply(countData_exp1_wo_outlier, as.numeric)
  class(countData_exp1_wo_outlier) <- "numeric"
  storage.mode(countData_exp1_wo_outlier) <- "numeric"}

# Import colData
colData <- read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Data/colData.tsv", header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
# Subset into Experiment groups:
colData_exp1_with_outlier <- subset(colData[c(1:7,22:28),])
colData_exp1_wo_outlier <- subset(colData_exp1_with_outlier[c(1:5, 7:9, 11:14),])

# design formula
designFormula <- "~ Clone + Group"

#  ---- With Outlier ----
## Build initial DEseq matrix
dds_exp1 <- DESeqDataSetFromMatrix(countData = countData_exp1_with_outlier, 
                                   colData = colData_exp1_with_outlier, 
                                   design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
dds_exp1 <- dds_exp1[ rowSums(DESeq2::counts(dds_exp1)) > 1, ]
# Now perform Differential expression analysis
dds_exp1 <- DESeq(dds_exp1)
# get results
DEresults_exp1_with_outlier <- as.data.frame(results(dds_exp1, contrast = c("Group", 'X', 'A')))


# ---- without outlier ----
## Build initial DEseq matrix
dds_exp1 <- DESeqDataSetFromMatrix(countData = countData_exp1_wo_outlier, 
                                   colData = colData_exp1_wo_outlier, 
                                   design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
dds_exp1 <- dds_exp1[ rowSums(DESeq2::counts(dds_exp1)) > 1, ]
# Now perform Differential expression analysis
dds_exp1 <- DESeq(dds_exp1)
# get results
DEresults_exp1_wo_outlier <- as.data.frame(results(dds_exp1, contrast = c("Group", 'X', 'A')))

# ---- intersect DEGs ----
sum(DEresults_exp1_with_outlier$padj < 0.1 & 
      (DEresults_exp1_with_outlier$log2FoldChange < 0 | DEresults_exp1_with_outlier$log2FoldChange > 0),
    na.rm = TRUE)
# 110 DEGs with outlier

sum(DEresults_exp1_wo_outlier$padj < 0.1 & 
      (DEresults_exp1_wo_outlier$log2FoldChange < 0 | DEresults_exp1_wo_outlier$log2FoldChange > 0),
    na.rm = TRUE)
# 211 WO outlier, as expected...

# extract list of DEGs:
str(DEresults_exp1_with_outlier)
str(DEresults_exp1_wo_outlier)

deg_with_outlier <- rownames(DEresults_exp1_with_outlier[
  !is.na(DEresults_exp1_with_outlier$padj) &
    DEresults_exp1_with_outlier$padj < 0.1 &
    DEresults_exp1_with_outlier$log2FoldChange != 0, 
])

deg_wo_outlier <- rownames(DEresults_exp1_wo_outlier[
  !is.na(DEresults_exp1_wo_outlier$padj) &
    DEresults_exp1_wo_outlier$padj < 0.1 &
    DEresults_exp1_wo_outlier$log2FoldChange != 0, 
])

length(deg_with_outlier)   # should be 110
length(deg_wo_outlier)     # should be 211

# Compare list of DEGs:
common_degs <- intersect(deg_with_outlier, deg_wo_outlier)
only_with_outlier <- setdiff(deg_with_outlier, deg_wo_outlier)
only_wo_outlier <- setdiff(deg_wo_outlier, deg_with_outlier)

length(common_degs)
length(only_with_outlier)
length(only_wo_outlier)

# optional: inspect them
head(common_degs)
head(only_with_outlier)
head(only_wo_outlier)

# Venn Diagram plot:
library(VennDiagram)
library(grid)

venn.plot <- venn.diagram(
  x = list(
    `With outlier` = deg_with_outlier,
    `Without outlier` = deg_wo_outlier
  ),
  filename = NULL,
  fill = c("steelblue", "tomato"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20)
)


png("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Results/Figures/Venn_of_Acute_with_and_without_outliers.png", width = 7, height = 5, units = "in", res = 600)

grid.newpage()
grid.draw(venn.plot)

dev.off()
