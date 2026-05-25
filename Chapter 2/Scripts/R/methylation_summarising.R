# Load required Packages
library(plyranges) # for data manipulation with genomic ranges  
library(GenomicRanges) # for converting and using files as GRanges  
library(qqman) # for Q-Q and Manhattan plots  
library(dplyr) # for data manipulation  
library(ggplot2) # for plotting  
library(data.table) # for data table manipulation
library(genomation) # for using annotateWithFeatures

# load data
## Sample Data
raw_cov <- fread("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/methylation_summarisation/meth_CpG_cov_reads")
colnames(raw_cov) <- c("chrom", "start", "end", "pc_meth", "count_meth", "count_un")
## Genome Data
genome <- read_gff3("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/methylation_summarisation/combined_annotations.gff3")

# Manipulate data to perform binomial analyses on methylation
cov_data <- raw_cov %>%
  mutate(
    coverage = count_meth + count_un
  ) %>%
  rowwise() %>%
  mutate(
    pv_meth = binom.test(
      count_meth, coverage, 0.01, alternative = "greater"
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    pv_meth_c = p.adjust(pv_meth, method = "fdr"),
    is_meth = ifelse(pv_meth_c < 0.05, 1, 0)
  )

# calculate mean methylation for all data in sample
mean_meth_rate =  sum(cov_data$is_meth)  / nrow(raw_cov)
write.table(mean_meth_rate,
            "mean_cpg_meth.txt") #overall methylation rate of Cs

# Change to granges object
cov_gr <- GRanges(
  seqnames = cov_data$chrom,
  ranges = IRanges(start = cov_data$start, end = cov_data$end),
  strand = "*",
  is_meth = cov_data$is_meth
)

## Annotate data with features (of interest):
# Annotate with all possible features:
{
region_types <- unique(mcols(genome)$type)
featuresList <- lapply(region_types, function(rt) {
  genome[genome$type == rt]
})
names(featuresList) <- region_types
featuresList <- GRangesList(featuresList)
}

region_types

# Annotate with features of interest:
{
promoters <- genome[genome$type == "promoter"]
genes <- genome[genome$type == "gene"]
exons <- genome[genome$type == "exon"]
downstream_regions <- genome[genome$type == "downstream_region"]
dispersed_repeats <- genome[genome$type == "dispersed_repeat"]
}

# Combine into a GRangesList
featuresList <- GRangesList(
  gene = genes,
  exon = exons,
  promoter = promoters,
  dispersed_repeat = dispersed_repeats,
  downstream_region = downstream_regions
)

# Annotate data with features
annotation_result <- annotateWithFeatures(
  target = cov_gr,
  features = featuresList,
  intersect.chr = TRUE
)

# Extract the stats
annotation_stats <- getTargetAnnotationStats(annotation_result, percentage = TRUE)
print(annotation_stats)

# annotation_stats might be returned as a named vector or data frame
# Ensure it's in a data frame format suitable for ggplot2:
annotation_df <- as.data.frame(annotation_stats)
colnames(annotation_df) <- "percentage"
# Convert percentage to numeric if it's not already
annotation_df$percentage <- as.numeric(annotation_df$percentage)
annotation_df$feature_type <- row.names(annotation_df)

# Create a pie chart using ggplot2:
ggplot(annotation_df, aes(x = "", y = percentage, fill = feature_type)) +
  geom_col(width = 1, color = "black") +
  coord_polar("y", start = 0) +
  # Add labels showing the percentage for each slice:
  geom_text(aes(label = paste0(percentage, "%")), 
            position = position_stack(vjust = 0.5, reverse = F), 
            size = 4, color = "black") +
  labs(
    title = "Distribution of CpGs - Sample 3-3",
    x = NULL,
    y = NULL,
    fill = "Feature Type"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )




# Make a plot of the stats (pie chart)
plotTargetAnnotation(annotation_result, precedence = TRUE, main = "Distribution of CpGs - Sample 3-3")
plot_filename <- paste0("TargetAnnotationPlot.png")
png(plot_filename, width = 1200, height = 800, res = 150)
plotTargetAnnotation(annotation_result, precedence = TRUE, main = "Distribution of CpGs - Sample 3-3")
dev.off()

# Calculating mean methylation per feature: ----

# Generate object (blank)
feature_assignment <- rep(NA, length(cov_gr))

# Assign features to each CpG
for (f in names(featuresList)) {
  # Find overlaps between targets and this feature type
  ov <- findOverlaps(cov_gr, featuresList[[f]])
  
  # Query hits correspond to indices in cov_gr
  q_hits <- queryHits(ov)
  
  # Only assign feature 'f' to those CpGs that are not yet assigned
  unassigned <- is.na(feature_assignment[q_hits])
  feature_assignment[q_hits[unassigned]] <- f
}

# Add the assigned feature type to cpgs_raw
raw_cov$feature_type <- feature_assignment

# Calculate mean methylation by feature type
mean_meth_rate_feature <- raw_cov %>%
  mutate(coverage = count_meth + count_un) %>%
  filter(coverage %in% 5:29) %>%
  group_by(feature_type) %>%
  summarise(
    n = n(),
    total_meth = sum(count_meth),
    mean = total_meth / n
  )

# Save the results
write.table(mean_meth_rate_feature,
            "mean_cpg_meth_by_feature.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


## Calculating the proportion of methylation cytosines per feature type: ----

library(GenomicRanges)
library(dplyr)

# Assume these objects already exist:
# cov_gr: a GRanges object containing your CpG data (including is_meth or other metadata)
# featuresList: a GRangesList, where each element is a different feature type

# Initialize a character vector to hold the assigned feature for each CpG
feature_assignment <- rep(NA_character_, length(cov_gr))

# Loop through each feature type in the featuresList
for (f in names(featuresList)) {
  # Find overlaps between cov_gr and this feature type
  ov <- findOverlaps(cov_gr, featuresList[[f]])
  
  # queryHits(ov) gives the indices of cov_gr that overlap with the current feature
  q_hits <- queryHits(ov)
  
  # Only assign feature 'f' to those CpGs that don't have a feature assigned yet
  unassigned <- is.na(feature_assignment[q_hits])
  feature_assignment[q_hits[unassigned]] <- f
}

# Add the assigned feature type as a metadata column to cov_gr
mcols(cov_gr)$feature_type <- feature_assignment

# Convert cov_gr to a data frame
cov_df <- as.data.frame(cov_gr)

# Now cov_df has columns for seqnames, start, end, width, strand, and any metadata columns 
# including 'feature_type' and any other columns that were originally in cov_gr (e.g. is_meth).
head(cov_df)

## Assuming cov_df includes:
# - "feature_type": The assigned feature type for each CpG (or NA if none)
# - "is_meth": A column indicating methylation status (1 for methylated, 0 for not)
# If "is_meth" is not numeric but logical/boolean, convert it to numeric: is_meth = as.numeric(is_meth)

prop_meth_by_feature <- cov_df %>%
  filter(!is.na(feature_type)) %>%  # Consider only CpGs assigned to a feature
  group_by(feature_type) %>%
  summarise(
    total_cpgs = n(),
    meth_cpgs = sum(is_meth, na.rm = TRUE),
    proportion_meth = meth_cpgs / total_cpgs
  )

# Inspect the summary
print(prop_meth_by_feature)

# Plot the proportion of methylated CpGs by feature type
ggplot(prop_meth_by_feature, aes(x = feature_type, y = proportion_meth, fill = feature_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Proportion of Methylated Cytosines by Feature Type",
    x = "Feature Type",
    y = "Proportion of Methylated CpGs"
  ) +
  coord_flip() # Flip coordinates for easier reading if you have many feature types
