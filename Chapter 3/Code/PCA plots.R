# Since we need to show the before and after of outlier removal, I will need to #
# load in the data again and then repeat some of the calculation steps

# Load a package list ----
# build function to quickly check packages
check_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      print(paste(pkg, "is not installed, installing package..."))
      install.packages(pkg)
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
      print(paste(pkg, "is installed and loaded"))
    }
  }
}

# list to load
pkg_list<-c("tidyverse", 
            "DESeq2", # DEA package
            "tibble", 
            "stats", 
            "EDASeq",
            "ggplot2", # plotting
            "pheatmap", 
            "dplyr", 
            "tidyr", 
            "ggfortify", 
            "corrplot", 
            "gprofiler2", 
            "knitr", 
            "gProfileR", 
            "gage", 
            "RUVSeq", 
            "ggforce",
            "here",
            "patchwork",
            "plotly",
            "plot3D",
            "geometry",
            "ggVennDiagram",
            "ggrepel",
            "gridExtra")

check_packages(pkg_list)

# Load data ----
# COunt Data
counts <- as.matrix(read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/Processed/A_Equina_Counts_unstranded_copy.txt", 
                               header = T, 
                               sep = '\t'))

# Load colData
colData <- read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/MS_colData.txt", header = T, sep = '\t', 
                      stringsAsFactors = TRUE)

# Calculate parameters ----
## CPM
cpm <- {apply(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)}
counts_rownames <- as.vector(row.names(counts))
rownames(cpm) <- counts_rownames

## Genelengths
geneLengths <- as.vector(subset(counts, select = c(Length)))
geneLengths <- as.vector(as.numeric(geneLengths))

## RPKM
rpkm <- {apply(X = subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)),
               MARGIN = 2, 
               FUN = function(x) {
                 10^9 * as.numeric(x) / geneLengths / sum(as.numeric(x))
               })
}
rownames(rpkm) <- counts_rownames

## RPK
rpk <- {apply(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/(geneLengths/1000))}
rownames(rpk) <- counts_rownames

## TPM
tpm <- apply(rpk, 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)
rownames(tpm) <- counts_rownames

## Variance
V <- apply(tpm, 1, var) # Can make Heatmaps from these values ... 

## Select top 100 genes:
selectedGenes <- names(V[order(V, decreasing = T)][1:100])

# PCA plots ----
## Function
make_pca_plot <- function(n_genes) {
  # 1) pick top genes
  sel <- names(sort(V, decreasing = TRUE))[seq_len(n_genes)]
  
  # 2) log-transform and make samples = rows
  M <- t(log2(tpm[sel, , drop = FALSE] + 1))
  
  # 3) PCA
  pca <- prcomp(M)  # samples as rows, genes as columns
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  
  # 4) join scores with metadata
  df <- as.data.frame(pca$x[, 1:2]) |>
    rownames_to_column("sample") |>
    left_join(
      colData |> as.data.frame() |> rownames_to_column("sample"),
      by = "sample"
    ) |>
    mutate(
      Group = recode(Group,
                     "A" = "Control",
                     "B" = "Diesel",
                     "C" = "Salinity",
                     "D" = "Diesel + Salinity"
      ),
      Group = factor(Group, levels = c("Control", "Diesel", "Salinity", "Diesel + Salinity"))
    )
  # 5) plot
  ggplot(df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3) +
    geom_mark_hull(aes(fill = Group), concavity = 5, expand = 0, radius = 0, alpha = 0.2, show.legend = FALSE) +
    labs(
      title = paste0("PCA (top ", n_genes, " genes)"),
      x = sprintf("PC1 (%.2f%%)", var_exp[1]),
      y = sprintf("PC2 (%.2f%%)", var_exp[2]),
      caption = "Counts: log2(x + 1)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold")) #+
  #geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf, min.segment.length = 0)
}

## PCAs with outlier ----

# Run for each target size
sizes <- c(100, 200, 500, 1000)
pca_plots <- setNames(map(sizes, make_pca_plot), paste0("top_", sizes))

# plot PCAs
p <- (pca_plots$top_100 | pca_plots$top_200) /
  (pca_plots$top_500 | pca_plots$top_1000)

p

ggsave("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Figures/Supplementary Figures/MS_PCA_with_outlier.png", 
       p, width = 16, height = 16)

# Removing outlier and repeating PCA plots ----
counts <- as.matrix(read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Data/Processed/A_Equina_Counts_unstranded_copy.txt", 
                               header = T, 
                               sep = '\t'))
counts_rm_outlier <- subset(counts, select = -X2) # removing sample "2" outlier
counts_rm_outlier_rownames <- as.vector(row.names(counts_rm_outlier))

# calculate parameters
## CPM
cpm <- {apply(subset(counts_rm_outlier, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)}
rownames(cpm) <- counts_rm_outlier_rownames

## GeneLengths
geneLengths <- as.vector(subset(counts_rm_outlier, select = c(Length)))
geneLengths <- as.vector(as.numeric(geneLengths))

## RPKM
rpkm <- {apply(X = subset(counts_rm_outlier, select = c(-Chr,-Start,-End,-Strand,-Length)),
               MARGIN = 2, 
               FUN = function(x) {
                 10^9 * as.numeric(x) / geneLengths / sum(as.numeric(x))
               })}
rownames(rpkm) <- counts_rm_outlier_rownames

## RPK
rpk <- {apply(subset(counts_rm_outlier, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/(geneLengths/1000))}
rownames(rpk) <- counts_rm_outlier_rownames

## TPM
tpm <- apply(rpk, 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)
rownames(tpm) <- counts_rm_outlier_rownames

# Vairnace
V <- apply(tpm, 1, var)

## Make plots ----
pca_plots_wo_outlier <- setNames(map(sizes, make_pca_plot), paste0("top_", sizes))

p <- (pca_plots_wo_outlier$top_100 | pca_plots_wo_outlier$top_200) /
  (pca_plots_wo_outlier$top_500 | pca_plots_wo_outlier$top_1000)

p

ggsave("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Figures/Supplementary Figures/MS_PCA_removed_outlier.png", 
       p, 
       width = 16, 
       height = 16)


# Main Figure plot
# plot PCAs with and without outlier together:
p2 <- (pca_plots$top_100 | pca_plots_wo_outlier$top_100)
p2

ggsave("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Figures/Main Figures/MS_PCA_outlier_comparison.png", 
       p2, 
       width = 16, 
       height = 8)


# Rmemove fluff and save workspace
rm(colData,
   counts,
   counts_rm_outlier,
   counts_rm_outlier_rownames,
   cpm,
   counts_rownames)

save.image("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/RData/PCA_plots.RData") 
load("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/RData/PCA_plots.RData")  
