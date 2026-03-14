# Script to correlate the LFCs of the LRT and Wald methods to see if they're comparable

# packages ----
## Make Function to Install and Load Package list
check_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      print(paste(pkg, "is not installed, installing package..."))
      install.packages(pkg) & suppressMessages(library(pkg, character.only = TRUE))
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
      print(paste(pkg, "is installed and loaded"))
    }
  }
}

## Generate a dataframe with a list of packages required
pkg_list<-c("DESeq2",
            "ggplot2",
            "dplyr",
            "tidyr",
            "ggrepel",
            "patchwork",
            "forcats",
            "purrr",
            "ComplexUpset",
            "grid",
            "ComplexHeatmap",
            "ggVennDiagram",
            "patchwork")

## Run Function::
check_packages(pkg_list)

# Load Data ----
# COunt Data
counts <- as.matrix(read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/Processed/A_Equina_Counts_unstranded_copy.txt", header = T, sep = '\t'))

countData <- as.matrix(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length,-X2)))
{apply(countData, 2, as.numeric)
  sapply(countData, as.numeric)
  class(countData) <- "numeric"
  storage.mode(countData) <- "numeric"}

# Load colData
colData <- read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/MS_colData.txt", header = T, sep = '\t', 
                      stringsAsFactors = TRUE)

# remove outlier "X2" from colData
colData <- subset(colData[c(1:24,26:32),])

reorder_idx <- match(colnames(countData), rownames(colData))
if (anyNA(reorder_idx)) {
  stop("Name mismatch after match(); check for typos/duplicates in sample names.\n",
       "Problem samples: ",
       paste(colnames(countData)[is.na(reorder_idx)], collapse=", "))
}
colData <- colData[reorder_idx, , drop = FALSE]

## Make the factors explicit and set "N" as the reference:
colData$Diesel   <- factor(colData$Diesel,   levels = c("N","Y"))
colData$Salinity <- factor(colData$Salinity, levels = c("N","Y"))

## LRT models ----
{
## Diesel ~ H0
# Subset data to only include samples with diesel or control treatments:

  keep <- colData$Salinity == "N"
  
  counData_keep <- countData[, keep]
  colData_keep <- droplevels(colData[keep,])
  
  dds <- DESeqDataSetFromMatrix(countData = countData[, keep],
                                colData   = droplevels(colData[keep,]),
                                design    = ~ Diesel)
  
  # (Optional but harmless) prefilter very low counts:
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  
  # LRT
  dds_diesel_H0 <- DESeq(dds, test = "LRT",
                         reduced = ~ 1)
  
  # Shrinking LFC for later comparison
  # see ~/Code/LRT_and_Wald_LFC_Comparison_plots.R for plotting
  res_diesel_H0_ashr <- lfcShrink(
    dds  = dds_diesel_H0,
    coef = "Diesel_Y_vs_N",   # replace with exact name from resultsNames(dds_diesel_H0)
    type = "ashr"
  )
  
  res_diesel_H0 <- results(dds_diesel_H0)
  res_diesel_H0 <- res_diesel_H0[order(res_diesel_H0$padj),]
  


## Salinity ~ H0
# Subset data to only include samples with Salinity or control treatments:

  keep <- colData$Diesel == "N"
  
  counData_keep <- countData[, keep]
  colData_keep <- droplevels(colData[keep,])
  
  dds <- DESeqDataSetFromMatrix(countData = countData[, keep],
                                colData   = droplevels(colData[keep,]),
                                design    = ~ Salinity)
  
  # (Optional but harmless) prefilter very low counts:
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  
  # LRT
  dds_Salinity_H0 <- DESeq(dds, test = "LRT",
                           reduced = ~ 1)
  
  # Shrinking LFC for later comparison
  # see ~/Code/LRT_and_Wald_LFC_Comparison_plots.R for plotting
  # shrink the Salinity coefficient using ashr for later
  res_salinity_H0_ashr <- lfcShrink(
    dds  = dds_Salinity_H0,
    coef = "Salinity_Y_vs_N",
    type = "ashr"
  )
  
  res_Salinity_H0 <- results(dds_Salinity_H0)
  res_Salinity_H0 <- res_Salinity_H0[order(res_Salinity_H0$padj),]
}

## Wald models ----
{
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData   = colData,
                              design    = ~ Diesel + Salinity + Diesel:Salinity)
# prefilter very low counts:
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ] 
dds_wald <- DESeq(dds, test = "Wald")
resultsNames(dds_wald)
# You'll likely see names like:
#   "Diesel_Y_vs_N"
#   "Salinity_Y_vs_N"
#   "DieselY.SalinityY"   # the interaction coefficient

## Diesel ~ H0
res_Diesel_at_SalN <- lfcShrink(dds_wald, coef = "Diesel_Y_vs_N", type = "ashr")
res_Diesel_at_SalN <- as.data.frame(res_Diesel_at_SalN)
## UNsrhunk results for comparison
res_Diesel_at_SalN_no_shrink <- results(dds_wald, name = "Diesel_Y_vs_N")
res_Diesel_at_SalN_no_shrink <- as.data.frame(res_Diesel_at_SalN_no_shrink)

# Salinity ~ H0
res_Sal_at_DieselN <- lfcShrink(dds_wald, coef = "Salinity_Y_vs_N", type = "ashr")
res_Sal_at_DieselN <- as.data.frame(res_Sal_at_DieselN)
# unshrunk results for comparison
res_Sal_at_DieselN_no_shrink <- results(dds_wald, name = "Salinity_Y_vs_N")
res_Sal_at_DieselN_no_shrink <- as.data.frame(res_Sal_at_DieselN_no_shrink)
}

# Data maniputlation ----
{
  # Diesel
diesel_lfc_df_raw <- data.frame(
  gene = rownames(res_diesel_H0),
  lfc_LRT_diesel = res_diesel_H0$log2FoldChange,
  stringsAsFactors = FALSE
)

diesel_wald_df_raw <- data.frame(
  gene = rownames(res_Diesel_at_SalN_no_shrink),
  lfc_Wald_diesel = res_Diesel_at_SalN_no_shrink$log2FoldChange,
  stringsAsFactors = FALSE
)


diesel_lfc_df_shrunk <- data.frame(
  gene = rownames(res_diesel_H0_ashr),
  lfc_LRT_diesel = res_diesel_H0_ashr$log2FoldChange,
  stringsAsFactors = FALSE
)

diesel_wald_df_shrunk <- data.frame(
  gene = rownames(res_Diesel_at_SalN),
  lfc_Wald_diesel = res_Diesel_at_SalN$log2FoldChange,
  stringsAsFactors = FALSE
)

diesel_cor_df <- list(
  data.frame(
    gene = rownames(res_diesel_H0_ashr),
    lfc_LRT_diesel_shrunk = res_diesel_H0_ashr$log2FoldChange
  ),
  data.frame(
    gene = rownames(res_Diesel_at_SalN),
    lfc_Wald_diesel_shrunk = res_Diesel_at_SalN$log2FoldChange
  ),
  data.frame(
    gene = rownames(res_diesel_H0),
    lfc_LRT_diesel_raw = res_diesel_H0$log2FoldChange
  ),
  data.frame(
    gene = rownames(res_Diesel_at_SalN_no_shrink),
    lfc_Wald_diesel_raw = res_Diesel_at_SalN_no_shrink$log2FoldChange
  )
) %>%
  purrr::reduce(inner_join, by = "gene") %>%
  filter(
    !is.na(lfc_LRT_diesel_shrunk),
    !is.na(lfc_Wald_diesel_shrunk),
    is.finite(lfc_LRT_diesel_shrunk),
    is.finite(lfc_Wald_diesel_shrunk),
    lfc_LRT_diesel_shrunk != 0,
    lfc_Wald_diesel_shrunk != 0,
    !is.na(lfc_LRT_diesel_raw),
    !is.na(lfc_Wald_diesel_raw),
    is.finite(lfc_LRT_diesel_raw),
    is.finite(lfc_Wald_diesel_raw),
    lfc_LRT_diesel_raw != 0,
    lfc_Wald_diesel_raw != 0
  )

# Salinity
sal_lfc_df <- data.frame(
  gene = rownames(res_salinity_H0_ashr),
  lfc_LRT_salinity = res_salinity_H0_ashr$log2FoldChange,
  stringsAsFactors = FALSE
)

sal_wald_df <- data.frame(
  gene = rownames(res_Sal_at_DieselN),
  lfc_Wald_salinity = res_Sal_at_DieselN$log2FoldChange,
  stringsAsFactors = FALSE
)

sal_cor_df <- inner_join(sal_lfc_df, sal_wald_df, by = "gene") %>%
  filter(
    !is.na(lfc_LRT_salinity),
    !is.na(lfc_Wald_salinity),
    is.finite(lfc_LRT_salinity),
    is.finite(lfc_Wald_salinity)
  ) %>%
  filter(
    !is.na(lfc_LRT_salinity),
    !is.na(lfc_Wald_salinity),
    is.finite(lfc_LRT_salinity),
    is.finite(lfc_Wald_salinity),
    lfc_LRT_salinity != 0,
    lfc_Wald_salinity != 0
  )

salinity_lfc_df_raw <- data.frame(
  gene = rownames(res_Salinity_H0),
  lfc_LRT_salinity = res_Salinity_H0$log2FoldChange,
  stringsAsFactors = FALSE
)

salinity_wald_df_raw <- data.frame(
  gene = rownames(res_Sal_at_DieselN_no_shrink),
  lfc_Wald_salinity = res_Sal_at_DieselN_no_shrink$log2FoldChange,
  stringsAsFactors = FALSE
)


salinity_lfc_df_shrunk <- data.frame(
  gene = rownames(res_salinity_H0_ashr),
  lfc_LRT_salinity = res_salinity_H0_ashr$log2FoldChange,
  stringsAsFactors = FALSE
)

salinity_wald_df_shrunk <- data.frame(
  gene = rownames(res_Sal_at_DieselN),
  lfc_Wald_salinity = res_Sal_at_DieselN$log2FoldChange,
  stringsAsFactors = FALSE
)

salinity_cor_df <- list(
  data.frame(
    gene = rownames(res_salinity_H0_ashr),
    lfc_LRT_salinity_shrunk = res_salinity_H0_ashr$log2FoldChange
  ),
  data.frame(
    gene = rownames(res_Sal_at_DieselN),
    lfc_Wald_salinity_shrunk = res_Sal_at_DieselN$log2FoldChange
  ),
  data.frame(
    gene = rownames(res_Salinity_H0),
    lfc_LRT_salinity_raw = res_Salinity_H0$log2FoldChange
  ),
  data.frame(
    gene = rownames(res_Sal_at_DieselN_no_shrink),
    lfc_Wald_salinity_raw = res_Sal_at_DieselN_no_shrink$log2FoldChange
  )
) %>%
  purrr::reduce(inner_join, by = "gene") %>%
  filter(
    !is.na(lfc_LRT_salinity_shrunk),
    !is.na(lfc_Wald_salinity_shrunk),
    is.finite(lfc_LRT_salinity_shrunk),
    is.finite(lfc_Wald_salinity_shrunk),
    lfc_LRT_salinity_shrunk != 0,
    lfc_Wald_salinity_shrunk != 0,
    !is.na(lfc_LRT_salinity_raw),
    !is.na(lfc_Wald_salinity_raw),
    is.finite(lfc_LRT_salinity_raw),
    is.finite(lfc_Wald_salinity_raw),
    lfc_LRT_salinity_raw != 0,
    lfc_Wald_salinity_raw != 0
  )
}

# Correlation Tests ----
# Pearson and Spearman correlations

## Diesel shrunk
cor.test(diesel_cor_df$lfc_LRT_diesel_shrunk, diesel_cor_df$lfc_Wald_diesel_shrunk, method = "pearson")
cor.test(diesel_cor_df$lfc_LRT_diesel_shrunk, diesel_cor_df$lfc_Wald_diesel_shrunk, method = "spearman")

## Diesel Raw
cor.test(diesel_cor_df$lfc_LRT_diesel_raw, diesel_cor_df$lfc_Wald_diesel_raw, method = "pearson")
cor.test(diesel_cor_df$lfc_LRT_diesel_raw, diesel_cor_df$lfc_Wald_diesel_raw, method = "spearman")

## Salinity Shrunk
cor.test(salinity_cor_df$lfc_LRT_salinity_shrunk, salinity_cor_df$lfc_Wald_salinity_shrunk, method = "pearson")
cor.test(salinity_cor_df$lfc_LRT_salinity_shrunk, salinity_cor_df$lfc_Wald_salinity_shrunk, method = "spearman")

# Salinity Raw
cor.test(salinity_cor_df$lfc_LRT_salinity_raw, salinity_cor_df$lfc_Wald_salinity_raw, method = "pearson")
cor.test(salinity_cor_df$lfc_LRT_salinity_raw, salinity_cor_df$lfc_Wald_salinity_raw, method = "spearman")

## Plots ----
diesel_LRT_vs_Wald_LFC_shrunk <- ggplot(diesel_cor_df, aes(x = lfc_LRT_diesel_shrunk, y = lfc_Wald_diesel_shrunk)) +
  geom_point(alpha = 0.4, size = 1.5) +
  labs(
    x = "LRT log2FoldChange - shrunk",
    y = "Wald log2FoldChange - shrunk"
  ) +
  theme_bw()  + ylim(-2,2)

diesel_LRT_vs_Wald_LFC_raw <- ggplot(diesel_cor_df, aes(x = lfc_LRT_diesel_raw, y = lfc_Wald_diesel_raw)) +
  geom_point(alpha = 0.4, size = 1.5) +
  labs(
    x = "LRT log2FoldChange",
    y = "Wald log2FoldChange"
  ) +
  theme_bw()  

salinity_LRT_vs_Wald_LFC_shrunk <- ggplot(salinity_cor_df, aes(x = lfc_LRT_salinity_shrunk, y = lfc_Wald_salinity_shrunk)) +
  geom_point(alpha = 0.4, size = 1.5) +
  labs(
    x = "LRT log2FoldChange - shrunk",
    y = "Wald log2FoldChange - shrunk"
  ) +
  theme_bw()  + xlim(-25,25) + ylim(-25,25)

salinity_LRT_vs_Wald_LFC_raw <- ggplot(salinity_cor_df, aes(x = lfc_LRT_salinity_raw, y = lfc_Wald_salinity_raw)) +
  geom_point(alpha = 0.4, size = 1.5) +
  labs(
    x = "LRT log2FoldChange",
    y = "Wald log2FoldChange"
  ) +
  theme_bw() + xlim(-25,25) + ylim(-25,25)

all_lfc_plots <-  (diesel_LRT_vs_Wald_LFC_raw | diesel_LRT_vs_Wald_LFC_shrunk) /
  (salinity_LRT_vs_Wald_LFC_raw | salinity_LRT_vs_Wald_LFC_shrunk) + 
  plot_layout(axis_titles = "collect")

all_lfc_plots

# save plot
ggsave("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Figures/Main Figures/LRT_and_Wald_LRT_correlation.png",  
       all_lfc_plots,
       width = 16, 
       height = 12)
