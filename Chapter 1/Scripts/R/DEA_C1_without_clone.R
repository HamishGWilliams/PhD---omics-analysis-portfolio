designFormula <- "~ Group"

countData <- as.matrix(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)))
countData_exp1 <- as.matrix(subset(counts_exp1, select = c(-Chr,-Start,-End,-Strand,-Length)))
countData_exp2 <- as.matrix(subset(counts_exp2, select = c(-Chr,-Start,-End,-Strand,-Length)))

# Here we tell R to make sure that it reads our data as numeric, as I had issues with this previously
{apply(countData, 2, as.numeric)
  sapply(countData, as.numeric)
  class(countData) <- "numeric"
  storage.mode(countData) <- "numeric"}

{apply(countData_exp1, 2, as.numeric)
  sapply(countData_exp1, as.numeric)
  class(countData_exp1) <- "numeric"
  storage.mode(countData_exp1) <- "numeric"}

{apply(countData_exp2, 2, as.numeric)
  sapply(countData_exp2, as.numeric)
  class(countData_exp2) <- "numeric"
  storage.mode(countData_exp2) <- "numeric"}

# removing clone5 (columns "X6" and "X24")

# Remove the 'width' column
countData_exp1 <- as.matrix(subset(counts_exp1, select = c(-Chr,-Start,-End,-Strand,-Length)))
# Remove outliers in Count Data
countData_exp1 <- subset(countData_exp1[,c(1:5, 7:9, 11:14)])
{apply(countData_exp1, 2, as.numeric)
  sapply(countData_exp1, as.numeric)
  class(countData_exp1) <- "numeric"
  storage.mode(countData_exp1) <- "numeric"}
# Remove outliers in colData Data
colData_exp1 <- subset(colData[c(1:7,22:28),])
colData_exp1 <- subset(colData_exp1[c(1:5, 7:9, 11:14),])

## Build initial DEseq matrix
dds_exp1 <- DESeqDataSetFromMatrix(countData = countData_exp1, 
                                   colData = colData_exp1, 
                                   design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
dds_exp1 <- dds_exp1[ rowSums(DESeq2::counts(dds_exp1)) > 1, ]
# Now perform Differential expression analysis
dds_exp1 <- DESeq(dds_exp1)

# Samples are used as the control group.
DEresults_exp1 = results(dds_exp1, contrast = c("Group", 'X', 'A'))
# Sort results by increasing p-value
DEresults_exp1 <- DEresults_exp1[order(DEresults_exp1$pvalue),]

# Summary of DEA
summary(DEresults_exp1)
# Total number of genes with p < 0.1
sum(DEresults_exp1$padj < 0.1 & (DEresults_exp1$log2FoldChange < 0 | DEresults_exp1$log2FoldChange > 0), na.rm=TRUE)
# 33
sum(DEresults_exp1$padj < 0.1 & (DEresults_exp1$log2FoldChange < 0), na.rm=TRUE)
# 18 down
sum(DEresults_exp1$padj < 0.1 & (DEresults_exp1$log2FoldChange > 0), na.rm=TRUE)
# 15 up

# Check if lfcshrink() is suitable - checking for low proportion of counts
{
  dds <- estimateSizeFactors(dds_exp1)
  normalized_counts <- counts(dds, normalized=TRUE)
  average_counts <- rowMeans(normalized_counts)
  low_count_threshold <- 10
  low_count_genes <- average_counts < low_count_threshold
  proportion_low_count_genes <- sum(low_count_genes) / length(average_counts)
  cat("Proportion of low-count genes:", proportion_low_count_genes, "\n")
  if (proportion_low_count_genes > 0.2) {
    cat("Consider using lfcShrink() due to a high proportion of low-count genes.\n")
  } else {
    cat("lfcShrink() may not be necessary based on the proportion of low-count genes.\n")
  }
}

# Applying LFC Shrink:
{
  # Find design of DESeq
  resultsNames(dds)
  # run lfcShrink()
  res_shrunk <- lfcShrink(dds, coef = "Group_X_vs_A", type = "ashr")
  # save df of results for later
  res_shrunk_exp1 <- as.data.frame(res_shrunk)
  # rlog for PCA plots
  dds_exp1_rlog <- rlog(dds_exp1)
}

# Summary of DEA
summary(res_shrunk_exp1)
# Total number of genes with p < 0.1
sum(res_shrunk_exp1$padj < 0.1 & (res_shrunk_exp1$log2FoldChange < 0 | res_shrunk_exp1$log2FoldChange > 0), na.rm=TRUE)
# 33
sum(res_shrunk_exp1$padj < 0.1 & (res_shrunk_exp1$log2FoldChange < 0), na.rm=TRUE)
# 18 down
sum(res_shrunk_exp1$padj < 0.1 & (res_shrunk_exp1$log2FoldChange > 0), na.rm=TRUE)
# 15 up

res_shrunk_exp1
str(res_shrunk_exp1)


df_shrunk <- as.data.frame(res_shrunk)  # Convert the shrunken results to a dataframe
df_shrunk$gene <- rownames(df_shrunk)  # Ensure gene identifiers are in a column

# Merge the two dataframes based on the gene identifiers
# Highlight Differentially Expressed Genes with Updated Criteria
# Highlight Differentially Expressed Genes with Updated Criteria
df_shrunk$diffexpressed <- 'P-adj > 0.1'
df_shrunk$diffexpressed[df_shrunk$padj <= 0.1] <- 'P-adj <= 0.1'

# Ensure diffexpressed is a factor and set the levels so "SHRUNK" is last
df_shrunk$diffexpressed <- factor(df_shrunk$diffexpressed, levels = c("P-adj > 0.1", "P-adj <= 0.1"))

# Plot
exp1_volcano_plot_without_outlier_with_lfcshrink <- ggplot(
  data = df_shrunk,
  aes(
    x = log2FoldChange,
    y = -log10(pvalue),
    colour = diffexpressed
  )
) + 
  geom_point(size = 1.5, alpha = 0.6) +
  geom_vline(xintercept = 0, col = "grey50", linetype = "dashed", lwd = 1) + 
  scale_colour_manual(values = c(
    "P-adj > 0.1" = "grey80",
    "P-adj <= 0.1" = "black"
  )) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(0,20)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5)) +
  theme_classic() +
  theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    x = expression("log"[2] * " Fold Change"),
    y = expression("-log"[10] * " p-value")
  ) +
  ggtitle(NULL)

# Save the merged results to find which genes are significant:
print(exp1_volcano_plot_without_outlier_with_lfcshrink)

png(
  "C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Results/Figures/Acute_DEA_WO_Clone.png",
  width = 4,
  height = 3,
  units = "in",
  res = 600
)
exp1_volcano_plot_without_outlier_with_lfcshrink
dev.off()
