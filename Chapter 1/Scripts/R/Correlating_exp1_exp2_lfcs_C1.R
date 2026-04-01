# Making a separate script since this needs to tweak a lot of stuff to get the right data

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
counts_exp2<- subset(counts[,c(1:5,13:26)])
# Remove the 'width' column
countData_exp1 <- as.matrix(subset(counts_exp1, select = c(-Chr,-Start,-End,-Strand,-Length)))
countData_exp2 <- as.matrix(subset(counts_exp2, select = c(-Chr,-Start,-End,-Strand,-Length)))
# Remove outliers in Count Data
countData_exp1 <- as.matrix(subset(countData_exp1[,c(1:5, 7:9, 11:14)]))

# Here we tell R to make sure that it reads our data as numeric, as I had issues with this previously
{apply(countData_exp1, 2, as.numeric)
  sapply(countData_exp1, as.numeric)
  class(countData_exp1) <- "numeric"
  storage.mode(countData_exp1) <- "numeric"}

{apply(countData_exp2, 2, as.numeric)
  sapply(countData_exp2, as.numeric)
  class(countData_exp2) <- "numeric"
  storage.mode(countData_exp2) <- "numeric"}

# Import colData
colData <- read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Data/colData.tsv", header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
# Subset into Experiment groups:
# Remove outliers in colData Data 1
colData_exp1 <- subset(colData[c(1:7,22:28),])
colData_exp1 <- subset(colData_exp1[c(1:5, 7:9, 11:14),])
colData_exp2<- subset(colData[8:21,])

# design formula
designFormula <- "~ Clone + Group"

# Experiment 1 ----
## Build initial DEseq matrix
dds_exp1 <- DESeqDataSetFromMatrix(countData = countData_exp1, 
                                   colData = colData_exp1, 
                                   design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
dds_exp1 <- dds_exp1[ rowSums(DESeq2::counts(dds_exp1)) > 1, ]
# Now perform Differential expression analysis
dds_exp1 <- DESeq(dds_exp1)

# No Shrinking applied
DEresults_exp1_no_shrinking <- as.data.frame(results(dds_exp1, contrast = c("Group", 'X', 'A')))

# Apply Shrinkging
dds <- estimateSizeFactors(dds_exp1)
DEresults_exp1_with_shrinking <- as.data.frame(lfcShrink(dds, coef = "Group_X_vs_A", type = "ashr"))
# rename lfc column to reflect shrinking
names(DEresults_exp1_with_shrinking)[names(DEresults_exp1_with_shrinking) == "log2FoldChange"] <- "log2FoldChange_Shrunk"

# Merge dfs
merged_DE_results_exp1 <- merge(DEresults_exp1_no_shrinking,
                                DEresults_exp1_with_shrinking,
                                by = 0)

# Make a new df containing only the data we need from this merged dataset
  # give better column names too

str(merged_DE_results_exp1)

Acute_Experiment_Data <- setNames(
                         merged_DE_results_exp1[, c(
                           "Row.names",
                           "log2FoldChange",
                           "log2FoldChange_Shrunk",
                           "padj.x"
                         )],
                         c(
                           "gene_id",
                           "log2FC_exp1",
                           "log2FC_shrunk_exp1",
                           "padj_exp1"
                         )
                       )

# Experiment 2 ----
## Build initial DEseq matrix
dds_exp2 <- DESeqDataSetFromMatrix(countData = countData_exp2, 
                                   colData = colData_exp2, 
                                   design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
dds_exp2 <- dds_exp2[ rowSums(DESeq2::counts(dds_exp2)) > 1, ]
# Now perform Differential expression analysis
dds_exp2 <- DESeq(dds_exp2)

# No Shrinking applied
DEresults_exp2_no_shrinking <- as.data.frame(results(dds_exp2, contrast = c("Group", 'C', 'B')))

# Apply Shrinkging
dds <- estimateSizeFactors(dds_exp2)
DEresults_exp2_with_shrinking <- as.data.frame(lfcShrink(dds, coef = "Group_C_vs_B", type = "ashr"))
# rename lfc column to reflect shrinking
names(DEresults_exp2_with_shrinking)[names(DEresults_exp2_with_shrinking) == "log2FoldChange"] <- "log2FoldChange_Shrunk"

# Merge dfs
merged_DE_results_exp2 <- merge(DEresults_exp2_no_shrinking,
                                DEresults_exp2_with_shrinking,
                                by = 0)

# Make a new df containing only the data we need from this merged dataset
# give better column names too

str(merged_DE_results_exp2)

Primed_Experiment_Data <- setNames(
  merged_DE_results_exp2[, c(
    "Row.names",
    "log2FoldChange",
    "log2FoldChange_Shrunk",
    "padj.x"
  )],
  c(
    "gene_id",
    "log2FC_exp2",
    "log2FC_shrunk_exp2",
    "padj_exp2"
  )
)

# Merge results together ----
merged_df <- merge(
  Acute_Experiment_Data,
  Primed_Experiment_Data,
  by = "gene_id",
  all = TRUE
)

# Correlate raw and shrunk lfcs of both experiments
  # select rows only where data is contained for all genes

# non-srhunk data
cor_df <- merged_df[complete.cases(merged_df[, c("log2FC_exp1", "log2FC_exp2")]), ]

cor_test <- cor.test(
  cor_df$log2FC_exp1,
  cor_df$log2FC_exp2,
  method = "pearson"
)

cor_test
cor_test$estimate   # correlation coefficient
cor_test$p.value    # p-value
cor_test$conf.int   # confidence interval


plot1 <- plot(
  cor_df$log2FC_exp1,
  cor_df$log2FC_exp2,
  xlab = "log2FC Exp1",
  ylab = "log2FC Exp2",
  main = "Correlation of log2FC between experiments",
  pch = 16
)


# shrunk data
cor_df <- merged_df[complete.cases(merged_df[, c("log2FC_shrunk_exp1", "log2FC_shrunk_exp2")]), ]

cor_test <- cor.test(
  cor_df$log2FC_shrunk_exp1,
  cor_df$log2FC_shrunk_exp2,
  method = "pearson"
)

cor_test
cor_test$estimate   # correlation coefficient
cor_test$p.value    # p-value
cor_test$conf.int   # confidence interval

plot2 <- plot(
  cor_df$log2FC_shrunk_exp1,
  cor_df$log2FC_shrunk_exp2,
  xlab = "Shrunk log2FC Exp1",
  ylab = "Shrunk log2FC Exp2",
  main = "Correlation of shrunk log2FC between experiments",
  pch = 16
) 


library(ggplot2)
library(gridExtra)

# Subset complete cases
cor_df_nonshrunk <- merged_df[
  complete.cases(merged_df[, c("log2FC_exp1", "log2FC_exp2")]), 
]

cor_df_shrunk <- merged_df[
  complete.cases(merged_df[, c("log2FC_shrunk_exp1", "log2FC_shrunk_exp2")]), 
]

# Correlation tests
cor_nonshrunk <- cor.test(
  cor_df_nonshrunk$log2FC_exp1,
  cor_df_nonshrunk$log2FC_exp2,
  method = "pearson"
)

cor_shrunk <- cor.test(
  cor_df_shrunk$log2FC_shrunk_exp1,
  cor_df_shrunk$log2FC_shrunk_exp2,
  method = "pearson"
)

# One common symmetric limit for both plots
common_limit <- max(abs(c(
  cor_df_nonshrunk$log2FC_exp1,
  cor_df_nonshrunk$log2FC_exp2,
  cor_df_shrunk$log2FC_shrunk_exp1,
  cor_df_shrunk$log2FC_shrunk_exp2
)), na.rm = TRUE)

# Non-shrunk plot
p1 <- ggplot(cor_df_nonshrunk, aes(x = log2FC_exp1, y = log2FC_exp2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, colour = "red") +
  coord_fixed(
    xlim = c(-common_limit, common_limit),
    ylim = c(-common_limit, common_limit)
  ) +
  labs(
    title = "Non-shrunk log2FC",
    x = "log2FC Exp1",
    y = "log2FC Exp2",
    subtitle = paste0(
      "Pearson r = ", round(cor_nonshrunk$estimate, 3),
      ", p = ", signif(cor_nonshrunk$p.value, 3)
    )
  ) +
  theme_bw()

# Shrunk plot
p2 <- ggplot(cor_df_shrunk, aes(x = log2FC_shrunk_exp1, y = log2FC_shrunk_exp2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, colour = "blue") +
  coord_fixed(
    xlim = c(-common_limit, common_limit),
    ylim = c(-common_limit, common_limit)
  ) +
  labs(
    title = "Shrunk log2FC",
    x = "log2FC_shrunk Exp1",
    y = "log2FC_shrunk Exp2",
    subtitle = paste0(
      "Pearson r = ", round(cor_shrunk$estimate, 3),
      ", p = ", signif(cor_shrunk$p.value, 3)
    )
  ) +
  theme_bw()

# Arrange side by side
grid.arrange(p1, p2, ncol = 2)

common_theme <- theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

p1 <- p1 + common_theme
p2 <- p2 + common_theme

grid.arrange(p1, p2, ncol = 2)

png("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Results/Figures/lfc_correlation_plots.png", width = 10, height = 5, units = "in", res = 600)
grid.arrange(p1, p2, ncol = 2)
dev.off()
