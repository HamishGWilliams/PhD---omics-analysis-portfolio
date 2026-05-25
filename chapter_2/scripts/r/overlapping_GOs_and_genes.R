
# Venn overlap of methylation vs expression enrichment results
# Acute / Primed
# Keeps full custom 4-way Venn plotting code

# load packages ----
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(readxl)
library(tidyr)
library(VennDiagram)
library(grid)


# 1. Paths and options # -------------------------------
project_dir <- "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio"

results_dir <- file.path(project_dir, "chapter_1/results/enrichment/venn_outputs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

use_only_significant <- TRUE
padj_cutoff <- 0.1

set_order <- c(
  "Acute Methylation",
  "Acute Expression",
  "Primed Methylation",
  "Primed Expression"
)

# 2. Find files # -------------------------------

search_dirs <- c(project_dir, getwd())
search_dirs <- search_dirs[dir.exists(search_dirs)]

all_files <- search_dirs %>%
  map(~ list.files(
    path = .x,
    pattern = "\\.(csv|xlsx|xls|tsv|txt)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )) %>%
  unlist() %>%
  unique()

file_index <- tibble(
  path = all_files,
  file_name = basename(all_files),
  file_name_lower = tolower(basename(all_files))
)

find_file <- function(must_contain, prefer_ext = NULL, required = TRUE) {
  hits <- file_index
  
  for (term in tolower(must_contain)) {
    hits <- hits %>%
      filter(str_detect(file_name_lower, fixed(term)))
  }
  
  if (!is.null(prefer_ext)) {
    preferred <- hits %>%
      filter(str_detect(file_name_lower, fixed(tolower(prefer_ext))))
    
    if (nrow(preferred) > 0) {
      hits <- preferred
    }
  }
  
  hits <- hits %>%
    arrange(nchar(path))
  
  if (nrow(hits) == 0) {
    msg <- paste("No file found for:", paste(must_contain, collapse = ", "))
    if (required) warning(msg)
    return(NA_character_)
  }
  
  if (nrow(hits) > 1) {
    message(
      "Multiple matches for ",
      paste(must_contain, collapse = ", "),
      "\nUsing: ",
      hits$path[1]
    )
  }
  
  hits$path[1]
}

file_manifest <- tibble(
  experiment = c(
    "Acute",
    "Primed",
    "Acute",
    "Primed",
    "Acute",
    "Primed"
  ),
  data_type = c(
    "Methylation",
    "Methylation",
    "Expression",
    "Expression",
    "Expression",
    "Expression"
  ),
  analysis_type = c(
    "ORA",
    "ORA",
    "ORA",
    "ORA",
    "GSEA",
    "GSEA"
  ),
  path = c(
    find_file(c("supplementary_table_ora_exp1_all_results"), ".xlsx"),
    find_file(c("supplementary_table_ora_exp2_all_results"), ".xlsx"),
    find_file(c("ora_exp1", "nested", "semantic"), ".csv"),
    find_file(c("ora_exp2", "nested", "semantic"), ".csv"),
    find_file(c("gsea_exp1", "semantic"), ".xlsx", required = FALSE),
    find_file(c("gsea_exp2", "semantic"), ".xlsx", required = FALSE)
  )
) %>%
  filter(!is.na(path), file.exists(path)) %>%
  mutate(set_name = paste(experiment, data_type))

file_manifest

# 3. Load files from manifest # -------------------------------

read_manifest_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  message("Reading: ", basename(path))
  
  switch(
    ext,
    "csv" = readr::read_csv(path, show_col_types = FALSE),
    "tsv" = readr::read_tsv(path, show_col_types = FALSE),
    "txt" = readr::read_delim(path, delim = "\t", show_col_types = FALSE),
    "xlsx" = readxl::read_excel(path),
    "xls" = readxl::read_excel(path),
    stop("Unsupported file extension: ", ext, " in file: ", path)
  )
}

loaded_results <- file_manifest %>%
  mutate(
    dataset_id = paste(experiment, data_type, analysis_type, sep = "_") %>%
      str_replace_all("\\s+", "_"),
    data = map(path, read_manifest_file)
  )

loaded_results


results_list <- loaded_results$data
names(results_list) <- loaded_results$dataset_id

results_list

acute_methylation_ora  <- results_list[["Acute_Methylation_ORA"]]
primed_methylation_ora <- results_list[["Primed_Methylation_ORA"]]

acute_expression_ora   <- results_list[["Acute_Expression_ORA"]]
primed_expression_ora  <- results_list[["Primed_Expression_ORA"]]

primed_expression_gsea <- results_list[["Primed_Expression_GSEA"]]

# 4. Subset to significant enrichment results # -------------------------------

get_padj_col <- function(df) {
  possible_cols <- c(
    "padj",
    "p.adjust",
    "p_adjust",
    "p.adj",
    "adjusted_p",
    "adjusted_p_value",
    "FDR",
    "qvalue",
    "q_value"
  )
  
  hit <- possible_cols[possible_cols %in% names(df)]
  
  if (length(hit) == 0) {
    stop(
      "No adjusted p-value column found. Available columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }
  
  hit[1]
}

filter_significant <- function(df, cutoff = padj_cutoff) {
  padj_col <- get_padj_col(df)
  
  df %>%
    filter(!is.na(.data[[padj_col]])) %>%
    filter(.data[[padj_col]] < cutoff)
}

acute_methylation_ora_sig <- filter_significant(acute_methylation_ora)
primed_methylation_ora_sig <- filter_significant(primed_methylation_ora)

acute_expression_ora_sig <- filter_significant(acute_expression_ora)
primed_expression_ora_sig <- filter_significant(primed_expression_ora)

primed_expression_gsea_sig <- filter_significant(primed_expression_gsea)

# check
sig_counts <- tibble(
  dataset = c(
    "Acute Methylation ORA",
    "Primed Methylation ORA",
    "Acute Expression ORA",
    "Primed Expression ORA",
    "Primed Expression GSEA"
  ),
  n_total = c(
    nrow(acute_methylation_ora),
    nrow(primed_methylation_ora),
    nrow(acute_expression_ora),
    nrow(primed_expression_ora),
    nrow(primed_expression_gsea)
  ),
  n_significant = c(
    nrow(acute_methylation_ora_sig),
    nrow(primed_methylation_ora_sig),
    nrow(acute_expression_ora_sig),
    nrow(primed_expression_ora_sig),
    nrow(primed_expression_gsea_sig)
  )
)

sig_counts

# 5. Remove duplicate rows and check uniqueness # -------------------------------

sig_dfs <- list(
  acute_methylation_ora_sig = acute_methylation_ora_sig,
  primed_methylation_ora_sig = primed_methylation_ora_sig,
  acute_expression_ora_sig = acute_expression_ora_sig,
  primed_expression_ora_sig = primed_expression_ora_sig,
  primed_expression_gsea_sig = primed_expression_gsea_sig
)

# Optional: check duplicate counts before removing
duplicate_check_before <- sig_dfs %>%
  imap_dfr(~ tibble(
    dataset = .y,
    n_rows = nrow(.x),
    n_unique_rows = nrow(distinct(.x)),
    n_duplicate_rows = nrow(.x) - nrow(distinct(.x)),
    all_rows_unique = nrow(.x) == nrow(distinct(.x))
  ))

duplicate_check_before

# all rows unique

# 6. Combine expression and methylation significant results # -------------------------------

add_result_origin <- function(df, experiment, method, data_type) {
  df %>%
    mutate(
      experiment = experiment,
      method = method,
      data_type = data_type,
      result_origin = paste(experiment, method, sep = "_"),
      .before = 1
    )
}

# Remove columns that are not needed before binding 

drop_unused_cols <- function(df) {
  df %>%
    select(-any_of(c("log2FoldChange")))
}


expression_dfs <- list(
  add_result_origin(
    acute_expression_ora_sig,
    experiment = "Acute",
    method = "ORA",
    data_type = "Expression"
  ),
  add_result_origin(
    primed_expression_ora_sig,
    experiment = "Primed",
    method = "ORA",
    data_type = "Expression"
  ),
  add_result_origin(
    primed_expression_gsea_sig,
    experiment = "Primed",
    method = "GSEA",
    data_type = "Expression"
  )
)

expression_dfs <- expression_dfs %>%
  map(drop_unused_cols)

combined_expression_sig <- bind_rows(expression_dfs) %>%
  distinct()


methylation_dfs <- list(
  add_result_origin(
    acute_methylation_ora_sig,
    experiment = "Acute",
    method = "ORA",
    data_type = "Methylation"
  ),
  add_result_origin(
    primed_methylation_ora_sig,
    experiment = "Primed",
    method = "ORA",
    data_type = "Methylation"
  )
)

methylation_dfs <- methylation_dfs %>%
  map(drop_unused_cols)

combined_methylation_sig <- bind_rows(methylation_dfs) %>%
  distinct()

combined_check <- tibble(
  dataset = c("combined_expression_sig", "combined_methylation_sig"),
  n_rows = c(
    nrow(combined_expression_sig),
    nrow(combined_methylation_sig)
  ),
  n_unique_rows = c(
    nrow(distinct(combined_expression_sig)),
    nrow(distinct(combined_methylation_sig))
  ),
  n_duplicate_rows = c(
    nrow(combined_expression_sig) - nrow(distinct(combined_expression_sig)),
    nrow(combined_methylation_sig) - nrow(distinct(combined_methylation_sig))
  ),
  all_rows_unique = c(
    nrow(combined_expression_sig) == nrow(distinct(combined_expression_sig)),
    nrow(combined_methylation_sig) == nrow(distinct(combined_methylation_sig))
  )
)

combined_check

# 7. Build 4 sets for Venn diagram # -------------------------------

stopifnot("Representative_GO" %in% names(combined_expression_sig))
stopifnot("semantic_cluster_description" %in% names(combined_methylation_sig))
stopifnot("experiment" %in% names(combined_expression_sig))
stopifnot("experiment" %in% names(combined_methylation_sig))

# 7. Build 4 sets for Venn diagram # -------------------------------

clean_terms <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    .[!is.na(.) & . != ""]
}

acute_expression_terms <- combined_expression_sig %>%
  filter(experiment == "Acute") %>%
  pull(`Representative_GO`) %>%
  clean_terms()

acute_methylation_terms <- combined_methylation_sig %>%
  filter(experiment == "Acute Experiment") %>%
  pull(semantic_cluster_description) %>%
  clean_terms()

primed_methylation_terms <- combined_methylation_sig %>%
  filter(experiment == "Primed Experiment") %>%
  pull(semantic_cluster_description) %>%
  clean_terms()

primed_expression_terms <- combined_expression_sig %>%
  filter(experiment == "Primed") %>%
  pull(`Representative_GO`) %>%
  clean_terms()

venn_list <- list(
  "Acute Expression" = acute_expression_terms,
  "Primed Expression" = primed_expression_terms,
  "Acute Methylation" = acute_methylation_terms,
  "Primed Methylation" = primed_methylation_terms


)

# Quick input check
venn_input_check <- tibble(
  set_name = names(venn_list),
  n_terms = purrr::map_int(venn_list, length)
)

venn_input_check

venn_plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#4C78A8", "#72B7B2", "#F58518", "#E45756"),
  alpha = 0.5,
  lwd = 2,
  col = "black",
  cex = 1.1,
  fontface = "bold",
  cat.cex = 1.1,
  cat.fontface = "bold",
  cat.dist = c(0.1, 0.1, 0.1, 0.1),
  margin = 0.01
)

grid.newpage()
grid.draw(venn_plot)


# 7a. Build raw Venn term lists without unique() # -------------------------------

clean_terms_raw <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    .[!is.na(.) & . != ""]
}

venn_list_raw <- list(
  "Acute Expression" = combined_expression_sig %>%
    filter(experiment == "Acute") %>%
    pull(`Representative_GO`) %>%
    clean_terms_raw(),
  
  "Acute Methylation" = combined_methylation_sig %>%
    filter(experiment == "Acute Experiment") %>%
    pull(semantic_cluster_description) %>%
    clean_terms_raw(),
  
  "Primed Methylation" = combined_methylation_sig %>%
    filter(experiment == "Primed Experiment") %>%
    pull(semantic_cluster_description) %>%
    clean_terms_raw(),
  
  "Primed Expression" = combined_expression_sig %>%
    filter(experiment == "Primed") %>%
    pull(`Representative_GO`) %>%
    clean_terms_raw()
)

# 7b. Identify duplicated GO / semantic terms within each set # -------------------------------

duplicated_venn_terms <- venn_list_raw %>%
  purrr::imap_dfr(~ tibble(
    set_name = .y,
    term = .x
  )) %>%
  count(set_name, term, name = "n_occurrences") %>%
  filter(n_occurrences > 1) %>%
  arrange(set_name, desc(n_occurrences), term)

duplicated_venn_terms

# 7c. Compare raw versus unique Venn input counts # -------------------------------

venn_input_check_raw_vs_unique <- venn_list_raw %>%
  purrr::imap_dfr(~ tibble(
    set_name = .y,
    n_raw_terms = length(.x),
    n_unique_terms = length(unique(.x)),
    n_duplicates_removed_by_unique = length(.x) - length(unique(.x))
  ))

venn_input_check_raw_vs_unique

write_csv(
  duplicated_venn_terms,
  file.path("duplicated_terms_removed_by_unique_in_venn_inputs.csv")
)


# 10. Build gene-level 4-way Venn inputs # -------------------------------

split_gene_column <- function(df, gene_col, sep_pattern, dataset_name) {
  if (!gene_col %in% names(df)) {
    stop("Column '", gene_col, "' not found in ", dataset_name)
  }
  
  df %>%
    select(all_of(gene_col)) %>%
    filter(!is.na(.data[[gene_col]])) %>%
    mutate(gene = as.character(.data[[gene_col]])) %>%
    separate_rows(gene, sep = sep_pattern) %>%
    mutate(gene = stringr::str_trim(gene)) %>%
    filter(!is.na(gene), gene != "") %>%
    pull(gene)
}

# Expression genes
acute_expression_genes_ora <- split_gene_column(
  acute_expression_ora_sig,
  gene_col = "geneID",
  sep_pattern = "/",
  dataset_name = "acute_expression_ora_sig"
)

primed_expression_genes_ora <- split_gene_column(
  primed_expression_ora_sig,
  gene_col = "geneID",
  sep_pattern = "/",
  dataset_name = "primed_expression_ora_sig"
)

primed_expression_genes_gsea <- split_gene_column(
  primed_expression_gsea_sig,
  gene_col = "core_enrichment",
  sep_pattern = "/",
  dataset_name = "primed_expression_gsea_sig"
)

# Methylation genes
acute_methylation_genes <- split_gene_column(
  acute_methylation_ora_sig,
  gene_col = "contributing_genes",
  sep_pattern = ";\\s*",
  dataset_name = "acute_methylation_ora_sig"
)

primed_methylation_genes <- split_gene_column(
  primed_methylation_ora_sig,
  gene_col = "contributing_genes",
  sep_pattern = ";\\s*",
  dataset_name = "primed_methylation_ora_sig"
)

acute_expression_genes <- acute_expression_genes_ora %>%
  unique()

primed_expression_genes <- c(
  primed_expression_genes_ora,
  primed_expression_genes_gsea
) %>%
  unique()

acute_methylation_genes <- acute_methylation_genes %>%
  unique()

primed_methylation_genes <- primed_methylation_genes %>%
  unique()

gene_venn_list <- list(
  "Acute Expression" = acute_expression_genes,
  "Primed Methylation" = primed_methylation_genes,
  "Acute Methylation" = acute_methylation_genes,
  "Primed Expression" = primed_expression_genes
)


gene_venn_input_check <- tibble(
  set_name = names(gene_venn_list),
  n_genes = purrr::map_int(gene_venn_list, length),
  n_unique_genes = purrr::map_int(gene_venn_list, ~ length(unique(.x))),
  all_genes_unique = n_genes == n_unique_genes
)

gene_venn_input_check

expression_gene_source_check <- tibble(
  source = c(
    "Acute Expression ORA",
    "Primed Expression ORA",
    "Primed Expression GSEA"
  ),
  n_raw_genes = c(
    length(acute_expression_genes_ora),
    length(primed_expression_genes_ora),
    length(primed_expression_genes_gsea)
  ),
  n_unique_genes = c(
    length(unique(acute_expression_genes_ora)),
    length(unique(primed_expression_genes_ora)),
    length(unique(primed_expression_genes_gsea))
  )
)

expression_gene_source_check

# 11. Plot gene-level 4-way Venn diagram # -------------------------------

gene_venn_plot <- venn.diagram(
  x = gene_venn_list,
  filename = NULL,
  fill = c("#4C78A8", "#72B7B2", "#F58518", "#E45756"),
  alpha = 0.5,
  lwd = 2,
  col = "black",
  cex = 1.1,
  fontface = "bold",
  cat.cex = 1.1,
  cat.fontface = "bold",
  cat.dist = c(0.1, 0.1, 0.1, 0.1),
  margin = 0.01
)

grid.newpage()
grid.draw(gene_venn_plot)



#### Overlapping signifcant DEGs with DMRs -----------------
# Gene overlap input construction
# Acute / Primed x Expression / Methylation

# Load packages ----
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(tidyr)
library(purrr)
library(tibble)

# 1. Paths and options ----

project_dir <- "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio"

expression_path <- file.path(
  project_dir,
  "chapter_1/results/model_outputs/DE_results_with_gene_namess.xlsx"
)

acute_methylation_path <- file.path(
  project_dir,
  "Chapter 2/Data/exp1_DMRs_gene_id_matched.txt"
)

primed_methylation_path <- file.path(
  project_dir,
  "Chapter 2/Data/exp2_DMRs_gene_id_matched.txt"
)

padj_cutoff <- 0.1
p_fdr_cutoff <- 0.1

# 2. Load files ----

expression_results <- readxl::read_excel(expression_path)

acute_methylation_results <- readr::read_tsv(
  acute_methylation_path,
  show_col_types = FALSE
)

primed_methylation_results <- readr::read_tsv(
  primed_methylation_path,
  show_col_types = FALSE
)

# 3. Helper functions ----

find_col <- function(df, candidates, df_name = "dataframe") {
  exact_hit <- candidates[candidates %in% names(df)]
  
  if (length(exact_hit) > 0) {
    return(exact_hit[1])
  }
  
  normalised_names <- names(df) %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_$", "")
  
  normalised_candidates <- candidates %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_$", "")
  
  idx <- match(normalised_candidates, normalised_names)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) {
    stop(
      "Could not find any of these columns in ", df_name, ": ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }
  
  names(df)[idx[1]]
}

normalise_experiment <- function(x) {
  x_clean <- x %>%
    as.character() %>%
    str_to_lower() %>%
    str_trim()
  
  case_when(
    str_detect(x_clean, "acute|exp1|experiment 1") ~ "Acute",
    str_detect(x_clean, "primed|exp2|experiment 2") ~ "Primed",
    TRUE ~ NA_character_
  )
}

clean_gene_vector <- function(x, sep_pattern = ";\\s*|/|,\\s*") {
  tibble(gene = as.character(x)) %>%
    filter(!is.na(gene)) %>%
    separate_rows(gene, sep = sep_pattern) %>%
    mutate(gene = str_trim(gene)) %>%
    filter(
      !is.na(gene),
      gene != "",
      !gene %in% c("NA", "NaN", "NULL", "character(0)")
    ) %>%
    distinct(gene) %>%
    pull(gene)
}

extract_expression_genes <- function(df, experiment_name, padj_cutoff = 0.1) {
  experiment_col <- find_col(
    df,
    candidates = c("Experiment", "experiment"),
    df_name = "expression_results"
  )
  
  padj_col <- find_col(
    df,
    candidates = c("padj", "p.adjust", "p_adjust", "p.adj", "FDR", "qvalue"),
    df_name = "expression_results"
  )
  
  gene_col <- find_col(
    df,
    candidates = c("Gene ID", "Gene_ID", "gene_id", "geneID", "gene"),
    df_name = "expression_results"
  )
  
  df %>%
    mutate(
      experiment_clean = normalise_experiment(.data[[experiment_col]]),
      padj_numeric = as.numeric(.data[[padj_col]])
    ) %>%
    filter(
      experiment_clean == experiment_name,
      !is.na(padj_numeric),
      padj_numeric < padj_cutoff
    ) %>%
    pull(all_of(gene_col)) %>%
    clean_gene_vector()
}

extract_methylation_genes <- function(df, experiment_name, p_fdr_cutoff = 0.1) {
  p_fdr_col <- find_col(
    df,
    candidates = c("p_fdr", "p.fdr", "pfdr", "FDR", "qvalue"),
    df_name = paste0(experiment_name, "_methylation_results")
  )
  
  gene_col <- find_col(
    df,
    candidates = c("gene_id_labelling", "gene_id_labeling", "gene_id", "geneID", "gene"),
    df_name = paste0(experiment_name, "_methylation_results")
  )
  
  df %>%
    mutate(
      p_fdr_numeric = as.numeric(.data[[p_fdr_col]])
    ) %>%
    filter(
      !is.na(p_fdr_numeric),
      p_fdr_numeric < p_fdr_cutoff
    ) %>%
    pull(all_of(gene_col)) %>%
    clean_gene_vector()
}

# 4. Extract significant gene lists ----

acute_expression_genes <- extract_expression_genes(
  expression_results,
  experiment_name = "Acute",
  padj_cutoff = padj_cutoff
)

primed_expression_genes <- extract_expression_genes(
  expression_results,
  experiment_name = "Primed",
  padj_cutoff = padj_cutoff
)

acute_methylation_genes <- extract_methylation_genes(
  acute_methylation_results,
  experiment_name = "Acute",
  p_fdr_cutoff = p_fdr_cutoff
)

primed_methylation_genes <- extract_methylation_genes(
  primed_methylation_results,
  experiment_name = "Primed",
  p_fdr_cutoff = p_fdr_cutoff
)

# 5. Store lists in required structure ----

gene_overlap_lists <- list(
  "Acute Expression" = acute_expression_genes,
  "Acute Methylation" = acute_methylation_genes,
  "Primed Methylation" = primed_methylation_genes,
  "Primed Expression" = primed_expression_genes
)

gene_overlap_lists

# 6. Check list sizes and uniqueness ----

gene_overlap_check <- tibble(
  set_name = names(gene_overlap_lists),
  n_genes = map_int(gene_overlap_lists, length),
  n_unique_genes = map_int(gene_overlap_lists, ~ length(unique(.x))),
  all_genes_unique = n_genes == n_unique_genes
)

gene_overlap_check

# 8. Plot 4-way gene overlap Venn diagram ----

library(VennDiagram)
library(grid)

# Ensure the order is correct
gene_overlap_lists <- gene_overlap_lists[c(
  "Acute Expression",
  "Primed Methylation",
  "Acute Methylation",
  "Primed Expression"
)]

# Create Venn object
gene_overlap_venn <- venn.diagram(
  x = gene_overlap_lists,
  filename = NULL,
  
  # Styling
  fill = c("#4C78A8", "#72B7B2", "#F58518", "#E45756"),
  alpha = 0.5,
  col = "black",
  lwd = 2,
  
  # Region labels
  cex = 1.1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Category labels
  cat.cex = 1.1,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = "black",
  
  # Layout
  margin = 0.01,
  cat.dist = c(0.1, 0.1, 0.1, 0.1)
)

# Draw plot in RStudio plot window
grid.newpage()
grid.draw(gene_overlap_venn)




# plot with all figures --------------------------

library(VennDiagram)
library(grid)
library(ggplot2)
library(ggplotify)
library(patchwork)
library(tibble)

# 1. Set group order and colours ----

# Desired visual / legend order
legend_order <- c(
  "Acute Expression",
  "Acute Methylation",
  "Primed Methylation",
  "Primed Expression"
)

# VennDiagram internal drawing order:
# 1 = lower-left
# 2 = lower-right
# 3 = middle-left
# 4 = middle-right
venn_order <- c(
  "Acute Expression",    # area1: lower-left
  "Primed Expression",   # area2: lower-right
  "Acute Methylation",   # area3: middle-left
  "Primed Methylation"   # area4: middle-right
)

venn_cols <- c(
  "Acute Expression"   = "#4C78A8",
  "Primed Methylation" = "#F58518",
  "Acute Methylation"  = "#72B7B2",
  "Primed Expression"  = "#E45756"
)

# 2. Store all Venn inputs and titles together ----

venn_inputs <- list(
  "A. DEG and DMR gene overlap" = gene_overlap_lists[venn_order],
  "B. Semantic GO term overlap" = venn_list[venn_order],
  "C. Enrichment-associated gene overlap" = gene_venn_list[venn_order]
)

# 3. Function to make one clean Venn panel ----

make_venn_panel <- function(x, title) {
  
  venn_obj <- venn.diagram(
    x = x,
    filename = NULL,
    fill = unname(venn_cols[names(x)]),
    alpha = 0.5,
    col = "black",
    lwd = 2,
    
    # overlap numbers
    cex = 1.6,
    fontface = "bold",
    
    # remove circle/category labels
    cat.cex = 0,
    cat.col = rep("transparent", length(x)),
    
    # slightly larger Venn inside panel
    margin = 0.03
  )
  
  venn_grob <- grobTree(
    children = do.call(gList, venn_obj)
  )
  
  as.ggplot(venn_grob) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = -2)
      ),
      plot.margin = margin(0, 2, 0, 2),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# 4. Build the three Venn panels ----

venn_panels <- Map(
  make_venn_panel,
  x = venn_inputs,
  title = names(venn_inputs)
)


# 5. Build shared legend only ----
venn_legend <- ggplot(
  tibble(group = factor(legend_order, levels = legend_order)),
  aes(x = group, y = 1, fill = group)
) +
  geom_point(
    shape = 21,
    size = 8,
    colour = "black",
    alpha = 0
  ) +
  scale_fill_manual(
    values = venn_cols,
    breaks = legend_order
  ) +
  guides(
    fill = guide_legend(
      title = "Groups",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 21,
        size = 9,
        colour = "black",
        alpha = 0.6
      )
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.1, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.margin = margin(t = -6, b = 0),
    plot.margin = margin(t = -10, b = 0),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# 6. Combine into one figure ----

combined_venn_figure <- wrap_plots(venn_panels, ncol = 3) /
  venn_legend +
  plot_layout(heights = c(12, 0.9))

combined_venn_figure

# 7. Save figure ----

ggsave(
  filename = file.path("combined_three_venn_plots_clean.png"),
  plot = combined_venn_figure,
  width = 16,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path("combined_three_venn_plots_clean.pdf"),
  plot = combined_venn_figure,
  width = 18,
  height = 8,
  bg = "white"
)


# 8. Extract exact intersection members from each Venn plot ----

library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(stringr)

# Use the human-readable order, not the VennDiagram draw order
group_order <- c(
  "Acute Expression",
  "Acute Methylation",
  "Primed Methylation",
  "Primed Expression"
)

clean_set_list <- function(x, group_order) {
  x[group_order] %>%
    map(~ .x %>%
          as.character() %>%
          str_trim() %>%
          .[!is.na(.) & . != ""] %>%
          unique() %>%
          sort())
}

get_exact_intersection <- function(set_list, groups) {
  
  # Items present in all selected groups
  in_selected <- Reduce(intersect, set_list[groups])
  
  # Remove items present in any non-selected group
  other_groups <- setdiff(names(set_list), groups)
  
  if (length(other_groups) > 0) {
    in_others <- set_list[other_groups] %>%
      unlist(use.names = FALSE) %>%
      unique()
    
    in_selected <- setdiff(in_selected, in_others)
  }
  
  sort(unique(in_selected))
}

make_exact_intersection_table <- function(set_list, plot_name, group_order) {
  
  set_list <- clean_set_list(set_list, group_order)
  
  group_combinations <- map(
    2:length(group_order),
    ~ combn(group_order, .x, simplify = FALSE)
  ) %>%
    flatten()
  
  group_combinations %>%
    map_dfr(function(groups) {
      
      members <- get_exact_intersection(set_list, groups)
      
      if (length(members) == 0) {
        return(tibble(
          plot_name = plot_name,
          intersection = paste(groups, collapse = " + "),
          n_groups = length(groups),
          n_items = 0,
          item = NA_character_
        ))
      }
      
      tibble(
        plot_name = plot_name,
        intersection = paste(groups, collapse = " + "),
        n_groups = length(groups),
        n_items = length(members),
        item = members
      )
    })
}

gene_overlap_intersections <- make_exact_intersection_table(
  set_list = gene_overlap_lists,
  plot_name = "DEG and DMR gene overlap",
  group_order = group_order
)

go_term_intersections <- make_exact_intersection_table(
  set_list = venn_list,
  plot_name = "Semantic GO term overlap",
  group_order = group_order
)

enrichment_gene_intersections <- make_exact_intersection_table(
  set_list = gene_venn_list,
  plot_name = "Enrichment-associated gene overlap",
  group_order = group_order
)

gene_overlap_intersections
go_term_intersections
enrichment_gene_intersections
