## Extra - Adding protein names to DEA results

# packages: ----
check_packages_ClusterProfileR <- function(pkg_list) {
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

# Generate a dataframe with a list of packages required

pkg_list_ClusterProfileR<-c("tidyverse", # manipulations
                            "tibble", # tables
                            "stats", # stats
                            "ggplot2", # plottin
                            "dplyr", # manipulation
                            "tidyr", # manipulation
                            "ggfortify", # plotting
                            "knitr",# markdown
                            "clusterProfiler", # enrichment analysis
                            "enrichplot", # plotting
                            "data.table", # data formatting
                            "plotly", # plotting
                            "ggVennDiagram", # plotting
                            "VennDiagram", # plotting
                            "limma", # plotting
                            "ggridges", # plotting
                            "stringr") # manipulations

# Run Function::
check_packages_ClusterProfileR(pkg_list_ClusterProfileR)

# main code ----

# Load in DE results files
DE_results_1 <- read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Results/DEA_results_exp1_no_outiler_and_shrunk_values.csv",
                           header = T, 
                           sep = ',', 
                           stringsAsFactors = TRUE)

DE_results_2 <- read.table("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Results/DEA_results_exp2_shrunk_values.csv",
                           header = T, 
                           sep = ',', 
                           stringsAsFactors = TRUE)

# only significant results
DE_results_1 <- as.data.frame(na.omit(DE_results_1[DE_results_1$padj <= 0.1,]))
DE_results_2 <- as.data.frame(na.omit(DE_results_2[DE_results_2$padj <= 0.1,]))

## Perform naming transformation to make gene IDs matching ----
DE_results_1$X <- gsub("\\.TU", "", DE_results_1$X) # Removing "TU" part of names
DE_results_1$X <- gsub("path", "mrna", DE_results_1$X) # Replacing "path" -> "mrna"
DE_results_2$X <- gsub("\\.TU", "", DE_results_2$X) # Removing "TU" part of names
DE_results_2$X <- gsub("path", "mrna", DE_results_2$X) # Replacing "path" -> "mrna"
names(DE_results_1)[names(DE_results_1) == "X"] <- "geneID"
names(DE_results_2)[names(DE_results_2) == "X"] <- "geneID"

# Keep only the best protein hit per gene (lowest e-value in V11)
matching_file_best <- matching_file %>%
  transmute(
    geneID = as.character(V1),
    proteinID = as.character(V13),
    evalue = as.numeric(V11)
  ) %>%
  group_by(geneID) %>%
  slice_min(order_by = evalue, n = 1, with_ties = FALSE) %>%
  ungroup()

# Named lookup: geneID -> proteinID
protein_map <- setNames(matching_file_best$proteinID, matching_file_best$geneID)

# Version 1: keep ORA rows intact, add protein IDs in geneID order
DE_results_1_with_protein_names <- DE_results_1 %>%
  mutate(
    geneID = as.character(geneID),
    proteinID = sapply(strsplit(geneID, "/"), function(genes) {
      paste(protein_map[genes], collapse = "/")
    })
  )

DE_results_2_with_protein_names <- DE_results_2 %>%
  mutate(
    geneID = as.character(geneID),
    proteinID = sapply(strsplit(geneID, "/"), function(genes) {
      paste(protein_map[genes], collapse = "/")
    })
  )

# re-order columns
DE_results_1_with_protein_names <- DE_results_1_with_protein_names %>%
  dplyr::select(geneID, proteinID, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  arrange(padj, desc(abs(log2FoldChange)))

DE_results_2_with_protein_names <- DE_results_2_with_protein_names %>%
  dplyr::select(geneID, proteinID, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  arrange(padj, desc(abs(log2FoldChange)))

# combined into 1 df
DE_results_1_with_protein_names$Experiment <- "Acute"
DE_results_2_with_protein_names$Experiment <- "Primed"

DE_results_with_protein_names <- rbind(DE_results_1_with_protein_names, 
                                       DE_results_2_with_protein_names)

DE_results_with_protein_names <- DE_results_with_protein_names %>%
  dplyr::select(Experiment, geneID, proteinID, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  arrange(Experiment, padj)

# save files for supplementary materials:
write.csv(DE_results_with_protein_names, "C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 1/Results/DE_results_with_protein_namess.csv", row.names = FALSE)