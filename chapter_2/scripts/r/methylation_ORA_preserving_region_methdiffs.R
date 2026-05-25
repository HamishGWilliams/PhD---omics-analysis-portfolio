# Load DMR data
dmr_data <- read.table(
  "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/chapter_2/data/temp_storage/exp1_DMRs_gene_id_matched.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# -------------------------------
# User-defined columns
# -------------------------------
region_col <- "region"
gene_col <- "gene_id_labelling"
sig_col <- "p_fdr"
sig_cutoff <- 0.1

# Change this if your methylation-difference column has a different name
methdiff_col <- "meth.diff"

# Optional context column
context_col <- "context"

if (!region_col %in% names(dmr_data)) stop("region_col not found in dmr_data")
if (!gene_col %in% names(dmr_data)) stop("gene_col not found in dmr_data")
if (!sig_col %in% names(dmr_data)) stop("sig_col not found in dmr_data")
if (!methdiff_col %in% names(dmr_data)) {
  stop("methdiff_col not found in dmr_data. Check names(dmr_data) and update methdiff_col.")
}

# -------------------------------
# Create a unique DMR / region-level identifier
# -------------------------------
dmr_data <- dmr_data %>%
  mutate(
    dmr_row_index = row_number(),
    region_clean = trimws(as.character(.data[[region_col]])),
    gene_clean = trimws(as.character(.data[[gene_col]])),
    methDiff_value = suppressWarnings(as.numeric(.data[[methdiff_col]])),
    sig_value = suppressWarnings(as.numeric(.data[[sig_col]])),
    context_value = if (context_col %in% names(.)) {
      as.character(.data[[context_col]])
    } else {
      NA_character_
    },
    dmr_feature_id = paste(
      "exp1",
      region_clean,
      gene_clean,
      dmr_row_index,
      sep = "__"
    )
  )

# Helpers ----

# Set all region types to analyse
region_types <- c(
  "promoter",
  "three_prime_UTR",
  "exon",
  "gene",
  "five_prime_UTR",
  "downstream_region",
  "dispersed_repeat",
  "ncRNA_gene"
)

# Helper function to convert x/y strings to numeric
convertFractionToNumeric <- function(x) {
  sapply(as.character(x), function(val) {
    parts <- strsplit(val, "/")[[1]]
    if (length(parts) == 2) {
      as.numeric(parts[1]) / as.numeric(parts[2])
    } else {
      NA_real_
    }
  })
}

# Helper function to process and expand ORA output
process_ora_results <- function(ora_obj, feature_map, ontology_label) {
  
  if (is.null(ora_obj) || is.null(ora_obj@result) || nrow(ora_obj@result) == 0) {
    return(NULL)
  }
  
  as_tibble(ora_obj@result) %>%
    mutate(
      geneID_gene_level = as.character(geneID),
      geneID = strsplit(as.character(geneID), "/")
    ) %>%
    tidyr::unnest(geneID) %>%
    rename(gene_id_labelling = geneID) %>%
    left_join(
      feature_map,
      by = "gene_id_labelling",
      relationship = "many-to-many"
    ) %>%
    group_by(ID) %>%
    mutate(
      n_unique_genes_for_term = n_distinct(gene_id_labelling),
      n_dmr_features_for_term = n_distinct(dmr_feature_id)
    ) %>%
    ungroup() %>%
    mutate(
      GO = ontology_label
    )
}

# revised ORA loop ----
# Store all results safely in a list
ora_results_exp1_by_region <- list()

for (target_region in region_types) {
  
  message("Processing region: ", target_region)
  
  # Subset data to selected region
  region_data <- dmr_data %>%
    filter(
      region_clean == target_region,
      !is.na(gene_clean),
      gene_clean != ""
    )
  
  if (nrow(region_data) == 0) {
    message("  No rows found for region: ", target_region)
    next
  }
  
  # Region-specific background universe:
  # all unique genes represented in this region category
  background_genes <- region_data %>%
    distinct(gene_clean) %>%
    pull(gene_clean)
  
  # Filter for significant DMRs
  sig_data <- region_data %>%
    filter(!is.na(sig_value), sig_value <= sig_cutoff)
  
  if (nrow(sig_data) == 0) {
    message("  No significant DMRs for region: ", target_region)
    next
  }
  
  # Feature-level map: preserves every significant DMR linked to each gene
  feature_map <- sig_data %>%
    transmute(
      gene_id_labelling = gene_clean,
      dmr_feature_id = dmr_feature_id,
      region_feature = region_clean,
      context = context_value,
      methDiff = methDiff_value,
      p_fdr = sig_value
    ) %>%
    filter(
      !is.na(gene_id_labelling),
      gene_id_labelling != ""
    )
  
  # ORA should use unique genes, not duplicated DMR rows
  gene_ids <- feature_map %>%
    distinct(gene_id_labelling) %>%
    pull(gene_id_labelling)
  
  if (length(gene_ids) == 0) {
    message("  No significant gene IDs found for region: ", target_region)
    next
  }
  
  message("  Significant DMR rows: ", nrow(sig_data))
  message("  Unique significant genes used for ORA: ", length(gene_ids))
  message("  Background genes for this region: ", length(background_genes))
  
  # Run ORA: BP
  ORA_BP <- tryCatch(
    clusterProfiler::enricher(
      gene = gene_ids,
      universe = background_genes,
      TERM2GENE = GENES_BP,
      TERM2NAME = TERMs_BP,
      minGSSize = 1,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {
      message("  BP ORA failed for ", target_region, ": ", conditionMessage(e))
      NULL
    }
  )
  
  # Run ORA: MF
  ORA_MF <- tryCatch(
    clusterProfiler::enricher(
      gene = gene_ids,
      universe = background_genes,
      TERM2GENE = GENES_MF,
      TERM2NAME = TERMs_MF,
      minGSSize = 1,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {
      message("  MF ORA failed for ", target_region, ": ", conditionMessage(e))
      NULL
    }
  )
  
  result_list <- list()
  
  bp_expanded <- process_ora_results(
    ora_obj = ORA_BP,
    feature_map = feature_map,
    ontology_label = "Biological Processes"
  )
  
  if (!is.null(bp_expanded)) {
    result_list[["BP"]] <- bp_expanded
  }
  
  mf_expanded <- process_ora_results(
    ora_obj = ORA_MF,
    feature_map = feature_map,
    ontology_label = "Molecular Functions"
  )
  
  if (!is.null(mf_expanded)) {
    result_list[["MF"]] <- mf_expanded
  }
  
  if (length(result_list) == 0) {
    message("  No ORA results found for region: ", target_region)
    next
  }
  
  # Combine and annotate results
  all_results <- bind_rows(result_list) %>%
    mutate(
      region = target_region,
      experiment = "Acute Experiment",
      GeneRatio_raw = GeneRatio,
      BgRatio_raw = BgRatio,
      GeneRatio = convertFractionToNumeric(GeneRatio),
      BgRatio = convertFractionToNumeric(BgRatio)
    )
  
  # Save to list
  ora_results_exp1_by_region[[target_region]] <- all_results
  
  # Optional: also save to environment using original object naming style
  assign(paste0("ORA_results_", target_region, "_exp1"), all_results)
  
  # Save to CSV
  output_file <- paste0("ORA_results_", target_region, "_exp1.csv")
  write.csv(all_results, output_file, row.names = FALSE)
  
  message("  Analysis complete and results saved for region: ", target_region)
}


# combine all avialable region results  ----
ORA_results_exp1_combined <- bind_rows(
  ora_results_exp1_by_region,
  .id = "region_from_list"
)

ORA_results_exp1_combined_clean <- unique(ORA_results_exp1_combined)

# important checks ----

ORA_results_exp1_combined_clean %>%
  select(
    ID,
    Description,
    GO,
    region,
    gene_id_labelling,
    dmr_feature_id,
    methDiff,
    p_fdr,
    context,
    Count,
    n_unique_genes_for_term,
    n_dmr_features_for_term
  ) %>%
  head()


# make plot ----
library(dplyr)
library(ggplot2)
library(stringr)

# -------------------------------
# Make term-level plotting data
# -------------------------------

ORA_results_exp1_plot_df <- ORA_results_exp1_combined %>%
  distinct(
    ID,
    Description,
    GO,
    region,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    .keep_all = FALSE
  ) %>%
  filter(!is.na(GeneRatio), !is.na(Count)) %>%
  arrange(GeneRatio)

# Optional: keep only significant ORA terms
ORA_results_exp1_plot_df_sig <- ORA_results_exp1_plot_df %>%
  filter(!is.na(p.adjust), p.adjust <= 0.1)

# If no significant terms exist, fall back to all terms
if (nrow(ORA_results_exp1_plot_df_sig) > 0) {
  ORA_plot_data <- ORA_results_exp1_plot_df_sig
} else {
  ORA_plot_data <- ORA_results_exp1_plot_df
}

# -------------------------------
# Plotting data
# -------------------------------

ORA_plot_data <- ORA_results_exp1_plot_df %>%
  filter(!is.na(p.adjust)) %>%
  group_by(GO, region) %>%
  arrange(p.adjust, desc(GeneRatio), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

ORA_plot_data <- ORA_plot_data[order(ORA_plot_data$GeneRatio), ]

ORA_exp1_dotplot <- ggplot(
  data = ORA_plot_data,
  aes(
    y = reorder(Description, Count),
    x = GeneRatio,
    color = region
  )
) +
  geom_point(
    aes(size = Count),
    stroke = 4
  ) +
  labs(
    title = "ORA results: Experiment 1",
    y = "GO Description",
    x = "Gene Ratio"
  ) +
  theme_light() +
  scale_size_continuous(
    breaks = seq(
      min(ORA_plot_data$Count, na.rm = TRUE),
      max(ORA_plot_data$Count, na.rm = TRUE),
      by = 1
    )
  ) +
  theme(
    legend.spacing = unit(0.5, "cm"),
    legend.key.height = unit(3, "lines"),
    axis.text.y = element_text(size = 11)
  ) +
  guides(size = guide_legend(title = "Count")) +
  scale_y_discrete(
    labels = function(x) str_wrap(x, width = 28)
  ) +
  facet_wrap(~GO, scales = "free")

plot(ORA_exp1_dotplot)


# version with boxplot ----
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# -------------------------------
# Region colours matched to previous figure
# -------------------------------

region_cols <- c(
  "3' UTR" = "#F2B84B",
  "5' UTR" = "#56B4E9",
  "Downstream" = "#009E73",
  "Exon" = "#8FD400",
  "Gene body" = "#0072B2",
  "ncRNA" = "#D95F02",
  "Promoter" = "#CC79A7",
  "Repeat" = "#FF1F1F"
)

# -------------------------------
# 1. Term-level ORA table
# -------------------------------

ORA_results_exp1_plot_df <- ORA_results_exp1_combined %>%
  distinct(
    ID,
    Description,
    GO,
    region,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    .keep_all = FALSE
  ) %>%
  filter(!is.na(GeneRatio), !is.na(Count)) %>%
  arrange(GeneRatio)

# Optional: keep only significant ORA terms
ORA_results_exp1_plot_df_sig <- ORA_results_exp1_plot_df %>%
  filter(!is.na(p.adjust), p.adjust <= 0.1)

# If no significant terms exist, fall back to all terms
if (nrow(ORA_results_exp1_plot_df_sig) > 0) {
  ORA_plot_base <- ORA_results_exp1_plot_df_sig
} else {
  ORA_plot_base <- ORA_results_exp1_plot_df
}

# -------------------------------
# 2. Select GO terms to display
# -------------------------------
# This keeps the top 10 terms per ontology and region,
# but the later semi_join keeps all regions for those selected GO IDs.

selected_go_terms <- ORA_plot_base %>%
  filter(!is.na(p.adjust)) %>%
  group_by(GO, region) %>%
  arrange(p.adjust, desc(GeneRatio), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  distinct(ID, Description, GO)

# To display literally all GO terms instead, use this instead:
# selected_go_terms <- ORA_plot_base %>%
#   distinct(ID, Description, GO)

# -------------------------------
# 3. Pull methylation differences for selected GO terms
# -------------------------------

ORA_methdiff_dot_df_fixed <- ORA_results_exp1_combined %>%
  mutate(
    methDiff = suppressWarnings(as.numeric(methDiff))
  ) %>%
  filter(
    !is.na(methDiff),
    !is.na(ID),
    !is.na(Description),
    !is.na(GO),
    !is.na(region),
    !is.na(dmr_feature_id)
  ) %>%
  semi_join(
    selected_go_terms,
    by = c("ID", "Description", "GO")
  ) %>%
  distinct(
    ID,
    Description,
    GO,
    region,
    dmr_feature_id,
    methDiff,
    p.adjust,
    Count,
    .keep_all = TRUE
  ) %>%
  mutate(
    region_label = recode(
      region,
      "three_prime_UTR" = "3' UTR",
      "five_prime_UTR" = "5' UTR",
      "downstream_region" = "Downstream",
      "exon" = "Exon",
      "gene" = "Gene body",
      "ncRNA_gene" = "ncRNA",
      "promoter" = "Promoter",
      "dispersed_repeat" = "Repeat",
      .default = region
    ),
    region_label = factor(
      region_label,
      levels = c(
        "3' UTR",
        "5' UTR",
        "Downstream",
        "Exon",
        "Gene body",
        "ncRNA",
        "Promoter",
        "Repeat"
      )
    ),
    GO_label = recode(
      GO,
      "Biological Processes" = "Biological process",
      "Molecular Functions" = "Molecular function",
      .default = GO
    ),
    
    # Critical fix:
    # make each GO term unique by ontology + GO ID + description
    term_key = paste(GO_label, ID, Description, sep = "___")
  )

# -------------------------------
# 4. Create unique nested GO labels
# -------------------------------

term_order <- ORA_methdiff_dot_df_fixed %>%
  group_by(GO_label, ID, Description, term_key) %>%
  summarise(
    min_padj = min(p.adjust, na.rm = TRUE),
    median_methDiff = median(methDiff, na.rm = TRUE),
    n_regions = n_distinct(region_label),
    n_dmr_features = n_distinct(dmr_feature_id),
    .groups = "drop"
  ) %>%
  arrange(GO_label, min_padj, desc(n_regions), desc(n_dmr_features)) %>%
  mutate(
    term_label = paste0(
      str_wrap(Description, width = 35),
      "\n",
      ID
    )
  )

# Use a named label vector so each unique term_key gets the right label
term_labels <- setNames(term_order$term_label, term_order$term_key)

# Reverse levels so most significant terms appear near the top
term_levels <- rev(term_order$term_key)

ORA_methdiff_dot_df_fixed <- ORA_methdiff_dot_df_fixed %>%
  mutate(
    term_key = factor(term_key, levels = term_levels)
  )

# -------------------------------
# 5. Plot: dots only, all regions retained per selected GO term
# -------------------------------

ORA_exp1_methdiff_dotplot_fixed <- ggplot(
  ORA_methdiff_dot_df_fixed,
  aes(
    x = methDiff,
    y = term_key,
    color = region_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_point(
    position = position_jitter(
      width = 0,
      height = 0.20,
      seed = 101
    ),
    alpha = 0.8,
    size = 2.5
  ) +
  scale_color_manual(
    values = region_cols,
    na.translate = FALSE
  ) +
  scale_y_discrete(
    labels = term_labels
  ) +
  labs(
    x = "Methylation difference",
    y = "GO term",
    color = "Genomic region"
  ) +
  facet_wrap(
    ~GO_label,
    scales = "free_y"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 8)),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 8.5, colour = "black", lineheight = 0.9),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.height = unit(0.7, "cm"),
    
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(8, 12, 8, 8)
  )

plot(ORA_exp1_methdiff_dotplot_fixed)
ggsave("ORA_exp1_methdiff_dotplot_fixed.png", ORA_exp1_methdiff_dotplot_fixed, width = 24, height = 16, units= "cm", dpi = 600)

library(dplyr)
library(stringr)

# -------------------------------
# Region and ontology labels
# -------------------------------

ORA_results_exp1_supplementary <- ORA_results_exp1_combined %>%
  mutate(
    region_label = recode(
      region,
      "three_prime_UTR" = "3' UTR",
      "five_prime_UTR" = "5' UTR",
      "downstream_region" = "Downstream",
      "exon" = "Exon",
      "gene" = "Gene body",
      "ncRNA_gene" = "ncRNA",
      "promoter" = "Promoter",
      "dispersed_repeat" = "Repeat",
      .default = region
    ),
    ontology = recode(
      GO,
      "Biological Processes" = "Biological process",
      "Molecular Functions" = "Molecular function",
      .default = GO
    ),
    methDiff = suppressWarnings(as.numeric(methDiff))
  ) %>%
  
  # No p.adjust filtering here: keeps significant and non-significant results
  group_by(
    experiment,
    ontology,
    region,
    region_label,
    ID,
    Description
  ) %>%
  summarise(
    GeneRatio = first(GeneRatio),
    BgRatio = first(BgRatio),
    GeneRatio_raw = if ("GeneRatio_raw" %in% names(.)) first(GeneRatio_raw) else NA,
    BgRatio_raw = if ("BgRatio_raw" %in% names(.)) first(BgRatio_raw) else NA,
    pvalue = first(pvalue),
    p.adjust = first(p.adjust),
    qvalue = first(qvalue),
    Count = first(Count),
    
    n_unique_genes = n_distinct(gene_id_labelling),
    n_dmr_features = n_distinct(dmr_feature_id),
    
    contributing_genes = paste(sort(unique(gene_id_labelling)), collapse = "; "),
    contributing_dmr_features = paste(sort(unique(dmr_feature_id)), collapse = "; "),
    
    methDiff_median = median(methDiff, na.rm = TRUE),
    methDiff_mean = mean(methDiff, na.rm = TRUE),
    methDiff_min = min(methDiff, na.rm = TRUE),
    methDiff_max = max(methDiff, na.rm = TRUE),
    methDiff_values = paste(methDiff[!is.na(methDiff)], collapse = "; "),
    
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      !is.na(p.adjust) & p.adjust <= 0.001 ~ "***",
      !is.na(p.adjust) & p.adjust <= 0.01 ~ "**",
      !is.na(p.adjust) & p.adjust <= 0.05 ~ "*",
      !is.na(p.adjust) & p.adjust <= 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(
    ontology,
    region_label,
    p.adjust,
    pvalue
  )

ORA_results_exp1_supplementary_clean <- ORA_results_exp1_supplementary %>%
  mutate(
    GeneRatio = if ("GeneRatio_raw" %in% names(.)) GeneRatio_raw else GeneRatio,
    BgRatio = if ("BgRatio_raw" %in% names(.)) BgRatio_raw else BgRatio
  ) %>%
  select(
    -any_of(c(
      "region_label",
      "GeneRatio_raw",
      "BgRatio_raw",
      "contributing_dmr_features",
      "methDiff_median",
      "methDiff_mean",
      "methDiff_min",
      "methDiff_max"
    ))
  )

write.csv(
  ORA_results_exp1_supplementary_clean,
  file = "Supplementary_Table_ORA_exp1_all_results.csv",
  row.names = FALSE
)

# install.packages("writexl") # only needed once
# library(writexl)
writexl::write_xlsx(
  ORA_results_exp1_supplementary_clean,
  path = "Supplementary_Table_ORA_exp1_all_results.xlsx"
)

## adding semantic clustering to table for supplementary table ----
library(dplyr)
library(stringr)
library(writexl)

# -------------------------------
# Prepare semantic lookup for joining
# -------------------------------

semantic_lookup_for_supp <- go_semantic_lookup_ora %>%
  transmute(
    ID = str_extract(as.character(ID), "GO:[0-9]{7}"),
    ontology_code = ontology,
    semantic_cluster_description = semantic_cluster_Description
  ) %>%
  distinct()

# -------------------------------
# Region and ontology labels
# -------------------------------

ORA_results_exp1_supplementary <- ORA_results_exp1_combined %>%
  mutate(
    ID = str_extract(as.character(ID), "GO:[0-9]{7}"),
    region_label = recode(
      region,
      "three_prime_UTR" = "3' UTR",
      "five_prime_UTR" = "5' UTR",
      "downstream_region" = "Downstream",
      "exon" = "Exon",
      "gene" = "Gene body",
      "ncRNA_gene" = "ncRNA",
      "promoter" = "Promoter",
      "dispersed_repeat" = "Repeat",
      .default = region
    ),
    ontology = recode(
      GO,
      "Biological Processes" = "Biological process",
      "Molecular Functions" = "Molecular function",
      .default = GO
    ),
    ontology_code = case_when(
      str_detect(tolower(ontology), "biol|process|^bp$") ~ "BP",
      str_detect(tolower(ontology), "mole|function|^mf$") ~ "MF",
      TRUE ~ NA_character_
    ),
    methDiff = suppressWarnings(as.numeric(methDiff))
  ) %>%
  
  # No p.adjust filtering here: keeps significant and non-significant results
  group_by(
    experiment,
    ontology,
    ontology_code,
    region,
    region_label,
    ID,
    Description
  ) %>%
  summarise(
    GeneRatio = first(GeneRatio),
    BgRatio = first(BgRatio),
    GeneRatio_raw = if ("GeneRatio_raw" %in% names(.)) first(GeneRatio_raw) else NA,
    BgRatio_raw = if ("BgRatio_raw" %in% names(.)) first(BgRatio_raw) else NA,
    pvalue = first(pvalue),
    p.adjust = first(p.adjust),
    qvalue = first(qvalue),
    Count = first(Count),
    
    n_unique_genes = n_distinct(gene_id_labelling),
    n_dmr_features = n_distinct(dmr_feature_id),
    
    contributing_genes = paste(sort(unique(gene_id_labelling)), collapse = "; "),
    contributing_dmr_features = paste(sort(unique(dmr_feature_id)), collapse = "; "),
    
    methDiff_median = median(methDiff, na.rm = TRUE),
    methDiff_mean = mean(methDiff, na.rm = TRUE),
    methDiff_min = min(methDiff, na.rm = TRUE),
    methDiff_max = max(methDiff, na.rm = TRUE),
    methDiff_values = paste(methDiff[!is.na(methDiff)], collapse = "; "),
    
    .groups = "drop"
  ) %>%
  left_join(
    semantic_lookup_for_supp,
    by = c("ID", "ontology_code")
  ) %>%
  mutate(
    significance = case_when(
      !is.na(p.adjust) & p.adjust <= 0.001 ~ "***",
      !is.na(p.adjust) & p.adjust <= 0.01 ~ "**",
      !is.na(p.adjust) & p.adjust <= 0.05 ~ "*",
      !is.na(p.adjust) & p.adjust <= 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(
    ontology,
    region_label,
    p.adjust,
    pvalue
  )

# -------------------------------
# Clean final supplementary table
# -------------------------------

ORA_results_exp1_supplementary_clean <- ORA_results_exp1_supplementary %>%
  mutate(
    GeneRatio = if ("GeneRatio_raw" %in% names(.)) GeneRatio_raw else GeneRatio,
    BgRatio = if ("BgRatio_raw" %in% names(.)) BgRatio_raw else BgRatio
  ) %>%
  select(
    -any_of(c(
      "ontology_code",
      "region_label",
      "GeneRatio_raw",
      "BgRatio_raw",
      "contributing_dmr_features",
      "methDiff_median",
      "methDiff_mean",
      "methDiff_min",
      "methDiff_max"
    ))
  )

# -------------------------------
# Export
# -------------------------------

write.csv(
  ORA_results_exp1_supplementary_clean,
  file = "Supplementary_Table_ORA_exp1_all_results.csv",
  row.names = FALSE
)

writexl::write_xlsx(
  ORA_results_exp1_supplementary_clean,
  path = "Supplementary_Table_ORA_exp1_all_results.xlsx"
)
