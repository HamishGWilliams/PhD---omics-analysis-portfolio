getwd()

# Load necessary libraries
library(ggplot2)
library(readr)

# HISAT2 ----

# Read the CSV file
hisat2_data <- read.csv("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Data/hisat2_merged_with_experimental_data_quotes.csv")
for (col in names(hisat2_data)) {
  hisat2_data[[col]] <- gsub('"', '', hisat2_data[[col]])
}
hisat2_data$OverallAlignmentRate <- as.numeric(hisat2_data$OverallAlignmentRate)

desired_order <- c("X","A","B","C")

hisat2_data$Group <-  factor(hisat2_data$Group, levels = desired_order)

hisat2_data <- hisat2_data %>%
  mutate(
    NamedGroup = recode(
      Group,
      "X" = "Acute - Control",
      "A" = "Acute - Treatment",
      "B" = "Primed - Treatment",
      "C" = "Primed - Control"
    )
  )

summary(hisat2_data[,c(1,3,5,6,7)])
by(hisat2_data$OverallAlignmentRate, hisat2_data$NamedGroup, summary)

# Plotting
par(mfrow = c(1,1))
alignment_boxplot <- expression(
  boxplot(hisat2_data$OverallAlignmentRate~hisat2_data$NamedGroup,
        xlab = "Group", ylab = "Overall Alignment Percentage (%)",
        main = "Overall alignment % of HISAT2 alignment ~ Group"
        )
)
eval(alignment_boxplot)

# T.tests for differences between experient treatments

# Check distributions of alignments
by(hisat2_data$OverallAlignmentRate, hisat2_data$Group, shapiro.test)
par(mfrow = c(2,2))
by(hisat2_data$OverallAlignmentRate, hisat2_data$Group, hist)

X_A <- hisat2_data[hisat2_data$Group == c("X","A"),]
t.test(X_A$OverallAlignmentRate~X_A$Group)
# data:  X_A$OverallAlignmentRate by X_A$Group
# t = 0.49779, df = 2.4414, p-value = 0.6599

X_B <- hisat2_data[hisat2_data$Group == c("X","B"),]
t.test(X_B$OverallAlignmentRate~X_B$Group)
# data:  X_B$OverallAlignmentRate by X_B$Group
# t = 0.67891, df = 3.7595, p-value = 0.5367

X_C <- hisat2_data[hisat2_data$Group == c("X","C"),]
t.test(X_C$OverallAlignmentRate~X_C$Group)
# data:  X_C$OverallAlignmentRate by X_C$Group
# t = 1.254, df = 3.1564, p-value = 0.2947

A_B <- hisat2_data[hisat2_data$Group == c("A","B"),]
t.test(A_B$OverallAlignmentRate~A_B$Group)
# data:  A_B$OverallAlignmentRate by A_B$Group
# t = 0.6223, df = 5.9811, p-value = 0.5567

A_C <- hisat2_data[hisat2_data$Group == c("A","C"),]
t.test(A_C$OverallAlignmentRate~A_C$Group)
# data:  A_C$OverallAlignmentRate by A_C$Group
# t = 1.3482, df = 4.9902, p-value = 0.2356

B_C <- hisat2_data[hisat2_data$Group == c("B","C"),]
t.test(B_C$OverallAlignmentRate~B_C$Group)
# data:  B_C$OverallAlignmentRate by B_C$Group
# t = 0.46843, df = 2.5193, p-value = 0.677


# FeatureCounts ----

# Read the CSV file
featurecounts_data <- read.csv("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Data/featurecounts_alignment_summary.csv")

# Boxplot and analyses of featureCounts alignment percent by Group
desired_order <- c("X","A","B","C")
featurecounts_data$Group <-  factor(featurecounts_data$Group, levels = desired_order)

featurecounts_data <- featurecounts_data %>%
  mutate(
    NamedGroup = recode(
      Group,
      "X" = "Acute - Control",
      "A" = "Acute - Treatment",
      "B" = "Primed - Treatment",
      "C" = "Primed - Control"
    )
  )

par(mfrow = c(1,1))
alignment_boxplot <- expression(
  boxplot(hisat2_data$OverallAlignmentRate~hisat2_data$NamedGroup,
          xlab = "Group", ylab = "Overall Alignment Percentage (%)",
          main = "Overall alignment % of HISAT2 alignment ~ Group",
          cex.axis = 0.6
  ))
  
eval(alignment_boxplot)

par(mfrow = c(1,1))
feature_counts_boxplot <- expression(
  boxplot(
    featurecounts_data$SuccessfullyAssignedPercent ~ featurecounts_data$NamedGroup,
    xlab = "Groups",
    ylab = "Aligned Percentage (%)",
    main = "Alignment % of featureCounts ~ Group",
    cex.axis = 0.6
  ))
eval(feature_counts_boxplot)


par(mfrow = c(1,2))
eval(alignment_boxplot)
eval(feature_counts_boxplot)
par(mfrow = c(1,1))

png(
  "C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Results/Figures/Alignment_and_counts_boxplots.png",
  width = 3200,
  height = 1800,
  res = 300
)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 1) + 0.1)

eval(alignment_boxplot)
eval(feature_counts_boxplot)

par(mfrow = c(1, 1))
dev.off()

by(featurecounts_data$SuccessfullyAssignedPercent, featurecounts_data$Group, shapiro.test)
par(mfrow = c(2,2))
by(featurecounts_data$SuccessfullyAssignedPercent, featurecounts_data$Group, hist)
par(mfrow = c(1,1))

# T.tests

X_A <- featurecounts_data[featurecounts_data$Group == c("X","A"),]
t.test(X_A$SuccessfullyAssignedPercent~X_A$Group)
# data:  X_A$SuccessfullyAssignedPercent by X_A$Group
# t = -0.56025, df = 5.2149, p-value = 0.5985

X_B <- featurecounts_data[featurecounts_data$Group == c("X","B"),]
t.test(X_B$SuccessfullyAssignedPercent~X_B$Group)
# data:  X_B$SuccessfullyAssignedPercent by X_B$Group
# t = 0.87428, df = 5.831, p-value = 0.4165

X_C <- featurecounts_data[featurecounts_data$Group == c("X","C"),]
t.test(X_C$SuccessfullyAssignedPercent~X_C$Group)
# data:  X_C$SuccessfullyAssignedPercent by X_C$Group
# t = 0.82691, df = 7.8425, p-value = 0.43277

A_B <- featurecounts_data[featurecounts_data$Group == c("A","B"),]
t.test(A_B$SuccessfullyAssignedPercent~A_B$Group)
# data:  A_B$SuccessfullyAssignedPercent by A_B$Group
# t = 0.98679, df = 4.9912, p-value = 0.3691

A_C <- featurecounts_data[featurecounts_data$Group == c("A","C"),]
t.test(A_C$SuccessfullyAssignedPercent~A_C$Group)
# data:  A_C$SuccessfullyAssignedPercent by A_C$Group
# t = 0.90224, df = 7.5644, p-value = 0.3947

B_C <- featurecounts_data[featurecounts_data$Group == c("B","C"),]
t.test(B_C$SuccessfullyAssignedPercent~B_C$Group)
# data:  B_C$SuccessfullyAssignedPercent by B_C$Group
# t = -0.70206, df = 7.8607, p-value = 0.5029