# Code for boxplots with semanticlaly clustered names, 
## ORA boxplot, showing both overrepresented and non-overrepresented semantically clustered descriptions
## keeping the non-overrepresented descriptions in since they are still significant results, and help to build a picture of regulatory targets for methylation from stress reponses
## but keeping information about overrepresentation to understand their represenatation in context with the hwole genomic background.

# version with semantically clustered GO labels ----
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(readxl)

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
# 0. Read ORA table with semantic cluster labels
# -------------------------------
# Uploaded file contains:
# Description = original GO term name
# semantic_cluster_description = semantically clustered GO term description

ORA_results_exp1_combined <- read_excel(
  "Supplementary_Table_ORA_exp1_all_results.xlsx"
) %>%
  rename(
    GO = ontology
  ) %>%
  mutate(
    ID = as.character(ID),
    Description = as.character(Description),
    semantic_cluster_description = as.character(semantic_cluster_description),
    GO = as.character(GO),
    region = as.character(region),
    
    pvalue = as.numeric(pvalue),
    p.adjust = as.numeric(p.adjust),
    qvalue = as.numeric(qvalue),
    Count = as.numeric(Count),
    n_dmr_features = as.numeric(n_dmr_features),
    
    # Use clustered label where available; otherwise fall back to original GO name
    Description_clustered = coalesce(
      na_if(str_squish(semantic_cluster_description), ""),
      Description
    )
  )

# -------------------------------
# 1. Term-level ORA table
# -------------------------------

ORA_results_exp1_plot_df <- ORA_results_exp1_combined %>%
  distinct(
    ID,
    Description,
    Description_clustered,
    GO,
    region,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    n_dmr_features,
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
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  distinct(ID, Description, Description_clustered, GO)

# To display literally all GO terms instead, use this instead:
# selected_go_terms <- ORA_plot_base %>%
#   distinct(ID, Description, Description_clustered, GO)

# -------------------------------
# 3. Pull methylation differences for selected GO terms
# -------------------------------
# In the uploaded supplementary tables, methylation differences are stored as
# semicolon-separated values in methDiff_values, rather than as one row per DMR.

ORA_methdiff_dot_df_fixed <- ORA_results_exp1_combined %>%
  mutate(
    methDiff_values = str_split(as.character(methDiff_values), ";\\s*")
  ) %>%
  unnest(methDiff_values) %>%
  mutate(
    methDiff = suppressWarnings(as.numeric(methDiff_values))
  ) %>%
  filter(
    !is.na(methDiff),
    !is.na(ID),
    !is.na(Description),
    !is.na(Description_clustered),
    !is.na(GO),
    !is.na(region)
  ) %>%
  semi_join(
    selected_go_terms,
    by = c("ID", "Description", "GO")
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
      "Biological process" = "Biological process",
      "Molecular Functions" = "Molecular function",
      "Molecular function" = "Molecular function",
      .default = GO
    ),
    
    # Keep GO ID in the key so clustered labels that are shared by
    # multiple GO IDs remain distinguishable.
    term_key = paste(GO_label, ID, Description_clustered, sep = "___")
  )

# -------------------------------
# 4. Create unique nested GO labels
# -------------------------------

term_order <- ORA_methdiff_dot_df_fixed %>%
  group_by(GO_label, ID, Description, Description_clustered, term_key) %>%
  summarise(
    min_padj = min(p.adjust, na.rm = TRUE),
    median_methDiff = median(methDiff, na.rm = TRUE),
    n_regions = n_distinct(region_label),
    n_methDiff_values = n(),
    .groups = "drop"
  ) %>%
  arrange(GO_label, min_padj, desc(n_regions), desc(n_methDiff_values)) %>%
  mutate(
    term_label = paste0(
      str_wrap(Description_clustered, width = 35),
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
    y = "Semantically clustered GO term",
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

ggsave(
  "ORA_exp1_methdiff_dotplot_semantic_clusters.png",
  ORA_exp1_methdiff_dotplot_fixed,
  width = 24,
  height = 16,
  units = "cm",
  dpi = 600
)

# facte_grid with semantic names: ----
# semantic-clustered GO methylation dotplot with all terms ----

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(readxl)

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
# 0. Read ORA table with semantic cluster labels
# -------------------------------

ORA_results_exp1_combined <- read_excel(
  "Supplementary_Table_ORA_exp1_all_results.xlsx"
) %>%
  rename(
    GO = ontology
  ) %>%
  mutate(
    ID = as.character(ID),
    Description = as.character(Description),
    semantic_cluster_description = as.character(semantic_cluster_description),
    GO = as.character(GO),
    region = as.character(region),
    
    pvalue = suppressWarnings(as.numeric(pvalue)),
    p.adjust = suppressWarnings(as.numeric(p.adjust)),
    qvalue = suppressWarnings(as.numeric(qvalue)),
    Count = suppressWarnings(as.numeric(Count)),
    n_dmr_features = suppressWarnings(as.numeric(n_dmr_features)),
    
    # Use semantic clustered term name where available
    Description_clustered = coalesce(
      na_if(str_squish(semantic_cluster_description), ""),
      Description
    )
  )

# -------------------------------
# 1. Term-level ORA table
# -------------------------------

padj_cutoff <- 0.1

ORA_results_exp1_plot_df <- ORA_results_exp1_combined %>%
  distinct(
    ID,
    Description,
    Description_clustered,
    GO,
    region,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    n_dmr_features,
    .keep_all = FALSE
  ) %>%
  filter(
    !is.na(ID),
    !is.na(Description),
    !is.na(Description_clustered),
    !is.na(GO),
    !is.na(region),
    !is.na(p.adjust)
  ) %>%
  mutate(
    overrepresentation_status = if_else(
      p.adjust <= padj_cutoff,
      "Significantly overrepresented",
      "Not significantly overrepresented"
    ),
    overrepresentation_status = factor(
      overrepresentation_status,
      levels = c(
        "Significantly overrepresented",
        "Not significantly overrepresented"
      )
    )
  )

# -------------------------------
# 2. Select GO terms to display
# -------------------------------
# All GO terms are retained.

selected_go_terms <- ORA_results_exp1_plot_df %>%
  distinct(
    ID,
    Description,
    Description_clustered,
    GO
  )

# -------------------------------
# 3. Pull methylation differences for all GO terms
# -------------------------------
# The supplementary table stores methylation differences as
# semicolon-separated values in methDiff_values.

ORA_methdiff_dot_df_fixed <- ORA_results_exp1_combined %>%
  mutate(
    methDiff_values = str_split(as.character(methDiff_values), ";\\s*")
  ) %>%
  unnest(methDiff_values) %>%
  mutate(
    methDiff = suppressWarnings(as.numeric(methDiff_values))
  ) %>%
  filter(
    !is.na(methDiff),
    !is.na(ID),
    !is.na(Description),
    !is.na(Description_clustered),
    !is.na(GO),
    !is.na(region),
    !is.na(p.adjust)
  ) %>%
  semi_join(
    selected_go_terms,
    by = c("ID", "Description", "GO")
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
      "Biological process" = "Biological process",
      "Molecular Functions" = "Molecular function",
      "Molecular function" = "Molecular function",
      .default = GO
    ),
    
    overrepresentation_status = if_else(
      p.adjust <= padj_cutoff,
      "Significantly overrepresented",
      "Not significantly overrepresented"
    ),
    overrepresentation_status = factor(
      overrepresentation_status,
      levels = c(
        "Significantly overrepresented",
        "Not significantly overrepresented"
      )
    ),
    
    # Keep GO ID in the key so repeated semantic cluster names
    # still remain unique.
    term_key = paste(
      overrepresentation_status,
      GO_label,
      ID,
      Description_clustered,
      sep = "___"
    )
  )

# -------------------------------
# 4. Create unique nested GO labels
# -------------------------------

term_order <- ORA_methdiff_dot_df_fixed %>%
  group_by(
    overrepresentation_status,
    GO_label,
    ID,
    Description,
    Description_clustered,
    term_key
  ) %>%
  summarise(
    min_padj = min(p.adjust, na.rm = TRUE),
    median_methDiff = median(methDiff, na.rm = TRUE),
    n_regions = n_distinct(region_label),
    n_methDiff_values = n(),
    .groups = "drop"
  ) %>%
  arrange(
    overrepresentation_status,
    GO_label,
    min_padj,
    desc(n_regions),
    desc(n_methDiff_values)
  ) %>%
  mutate(
    term_label = paste0(
      str_wrap(Description_clustered, width = 35),
      "\n",
      ID
    )
  )

term_labels <- setNames(term_order$term_label, term_order$term_key)

term_levels <- rev(term_order$term_key)

ORA_methdiff_dot_df_fixed <- ORA_methdiff_dot_df_fixed %>%
  mutate(
    term_key = factor(term_key, levels = term_levels)
  )

# -------------------------------
# 5. Plot all terms
# -------------------------------

ORA_exp1_methdiff_dotplot_semantic_all_terms <- ggplot(
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
    labels = term_labels,
    drop = TRUE
  ) +
  labs(
    x = "Methylation difference",
    y = "Semantically clustered GO term",
    color = "Genomic region"
  ) +
  facet_wrap(
    overrepresentation_status~GO_label,
    scales = "free_y",
    space = "free_y"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 8)),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 7.5, colour = "black", lineheight = 0.85),
    
    strip.background = element_blank(),
    strip.text.y = element_text(size = 11, face = "bold", angle = 0),
    
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

plot(ORA_exp1_methdiff_dotplot_semantic_all_terms)

ggsave(
  "ORA_exp1_methdiff_dotplot_semantic_clusters_all_terms_sig_vs_nonsig.png",
  ORA_exp1_methdiff_dotplot_semantic_all_terms,
  width = 30,
  height = 40,
  units = "cm",
  dpi = 600
)

## horizontal boxplots instead ----
# semantic-clustered GO methylation horizontal boxplot with all terms ----

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(readxl)

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
# 0. Read ORA table with semantic cluster labels
# -------------------------------

ORA_results_exp1_combined <- read_excel(
  "Supplementary_Table_ORA_exp1_all_results.xlsx"
) %>%
  rename(
    GO = ontology
  ) %>%
  mutate(
    ID = as.character(ID),
    Description = as.character(Description),
    semantic_cluster_description = as.character(semantic_cluster_description),
    GO = as.character(GO),
    region = as.character(region),
    
    pvalue = suppressWarnings(as.numeric(pvalue)),
    p.adjust = suppressWarnings(as.numeric(p.adjust)),
    qvalue = suppressWarnings(as.numeric(qvalue)),
    Count = suppressWarnings(as.numeric(Count)),
    n_dmr_features = suppressWarnings(as.numeric(n_dmr_features)),
    
    # Use semantic clustered term name where available
    Description_clustered = coalesce(
      na_if(str_squish(semantic_cluster_description), ""),
      Description
    )
  )

# -------------------------------
# 1. Term-level ORA table
# -------------------------------

padj_cutoff <- 0.1

ORA_results_exp1_plot_df <- ORA_results_exp1_combined %>%
  distinct(
    ID,
    Description,
    Description_clustered,
    GO,
    region,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    n_dmr_features,
    .keep_all = FALSE
  ) %>%
  filter(
    !is.na(ID),
    !is.na(Description),
    !is.na(Description_clustered),
    !is.na(GO),
    !is.na(region),
    !is.na(p.adjust)
  ) %>%
  mutate(
    overrepresentation_status = if_else(
      p.adjust <= padj_cutoff,
      "Significantly overrepresented",
      "Not significantly overrepresented"
    ),
    overrepresentation_status = factor(
      overrepresentation_status,
      levels = c(
        "Significantly overrepresented",
        "Not significantly overrepresented"
      )
    )
  )

# -------------------------------
# 2. Select GO terms to display
# -------------------------------
# All GO terms are retained.

selected_go_terms <- ORA_results_exp1_plot_df %>%
  distinct(
    ID,
    Description,
    Description_clustered,
    GO
  )

# -------------------------------
# 3. Pull methylation differences for all GO terms
# -------------------------------
# The supplementary table stores methylation differences as
# semicolon-separated values in methDiff_values.

ORA_methdiff_boxplot_df <- ORA_results_exp1_combined %>%
  mutate(
    methDiff_values = str_split(as.character(methDiff_values), ";\\s*")
  ) %>%
  unnest(methDiff_values) %>%
  mutate(
    methDiff = suppressWarnings(as.numeric(methDiff_values))
  ) %>%
  filter(
    !is.na(methDiff),
    !is.na(ID),
    !is.na(Description),
    !is.na(Description_clustered),
    !is.na(GO),
    !is.na(region),
    !is.na(p.adjust)
  ) %>%
  semi_join(
    selected_go_terms,
    by = c("ID", "Description", "GO")
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
      "Biological process" = "Biological process",
      "Molecular Functions" = "Molecular function",
      "Molecular function" = "Molecular function",
      .default = GO
    ),
    
    overrepresentation_status = if_else(
      p.adjust <= padj_cutoff,
      "Significantly overrepresented",
      "Not significantly overrepresented"
    ),
    overrepresentation_status = factor(
      overrepresentation_status,
      levels = c(
        "Significantly overrepresented",
        "Not significantly overrepresented"
      )
    ),
    
    # Keep the hidden key unique by including GO ID,
    # but do not show the GO ID in the plotted axis label.
    term_key = paste(
      overrepresentation_status,
      GO_label,
      ID,
      Description_clustered,
      sep = "___"
    )
  )

# -------------------------------
# 4. Create unique nested semantic labels
# -------------------------------

term_order <- ORA_methdiff_boxplot_df %>%
  group_by(
    overrepresentation_status,
    GO_label,
    ID,
    Description,
    Description_clustered,
    term_key
  ) %>%
  summarise(
    min_padj = min(p.adjust, na.rm = TRUE),
    median_methDiff = median(methDiff, na.rm = TRUE),
    n_regions = n_distinct(region_label),
    n_methDiff_values = n(),
    .groups = "drop"
  ) %>%
  arrange(
    overrepresentation_status,
    GO_label,
    min_padj,
    desc(n_regions),
    desc(n_methDiff_values)
  ) %>%
  mutate(
    # Only semantic descriptions are shown on the plot
    term_label = str_wrap(Description_clustered, width = 35)
  )

term_labels <- setNames(term_order$term_label, term_order$term_key)

term_levels <- rev(term_order$term_key)

ORA_methdiff_boxplot_df <- ORA_methdiff_boxplot_df %>%
  mutate(
    term_key = factor(term_key, levels = term_levels)
  )

# -------------------------------
# 5. Horizontal boxplot
# -------------------------------

ORA_exp1_methdiff_boxplot_semantic_all_terms <- ggplot(
  ORA_methdiff_boxplot_df,
  aes(
    x = methDiff,
    y = term_key,
    fill = region_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_boxplot(
    orientation = "y",
    position = position_dodge2(
      width = 0.75,
      preserve = "single"
    ),
    width = 0.65,
    alpha = 0.8,
    outlier.size = 0.9,
    outlier.alpha = 0.9,
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = region_cols,
    na.translate = FALSE
  ) +
  scale_y_discrete(
    labels = term_labels,
    drop = TRUE
  ) +
  labs(
    x = "Methylation difference",
    y = "Semantically clustered GO term",
    fill = "Genomic region"
  ) +
  facet_wrap(overrepresentation_status~GO_label,
    scales = "free_y",
    space = "free_y"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 8)),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 7.5, colour = "black", lineheight = 0.85),
    
    strip.background = element_blank(),
    strip.text.y = element_text(size = 11, face = "bold", angle = 0),
    
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

plot(ORA_exp1_methdiff_boxplot_semantic_all_terms)

ggsave(
  "ORA_exp1_methdiff_boxplot_semantic_clusters_all_terms_sig_vs_nonsig.png",
  ORA_exp1_methdiff_boxplot_semantic_all_terms,
  width = 30,
  height = 40,
  units = "cm",
  dpi = 600
)

# better plot ----
# -------------------------------
# Collapse GO terms with the same semantic description into one y-axis row
# -------------------------------

term_order_semantic <- ORA_methdiff_boxplot_df %>%
  group_by(Description_clustered) %>%
  summarise(
    min_padj = min(p.adjust, na.rm = TRUE),
    median_methDiff = median(methDiff, na.rm = TRUE),
    n_regions = n_distinct(region_label),
    n_methDiff_values = n(),
    .groups = "drop"
  ) %>%
  arrange(min_padj, desc(n_regions), desc(n_methDiff_values)) %>%
  mutate(
    semantic_term_key = Description_clustered,
    semantic_term_label = str_wrap(Description_clustered, width = 35)
  )

semantic_term_labels <- setNames(
  term_order_semantic$semantic_term_label,
  term_order_semantic$semantic_term_key
)

semantic_term_levels <- rev(term_order_semantic$semantic_term_key)

ORA_methdiff_boxplot_df <- ORA_methdiff_boxplot_df %>%
  mutate(
    semantic_term_key = factor(
      Description_clustered,
      levels = semantic_term_levels
    )
  )

# -------------------------------
# Horizontal boxplot + points
# Facets are columns.
# Same semantic descriptions share one row.
# -------------------------------

ORA_exp1_methdiff_boxplot_semantic_all_terms <- ggplot(
  ORA_methdiff_boxplot_df,
  aes(
    x = methDiff,
    y = semantic_term_key,
    fill = region_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_boxplot(
    orientation = "y",
    position = position_dodge2(
      width = 0.75,
      preserve = "single"
    ),
    width = 0.65,
    alpha = 0.65,
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  geom_point(
    aes(
      colour = region_label
    ),
    position = position_jitterdodge(
      dodge.width = 0.75,
      jitter.width = 0,
      jitter.height = 0.08,
      seed = 101
    ),
    size = 0.9,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = region_cols,
    na.translate = FALSE
  ) +
  scale_colour_manual(
    values = region_cols,
    na.translate = FALSE
  ) +
  #scale_y_discrete(
  #  labels = semantic_term_labels,
  #  drop = FALSE
  #) +
  labs(
    x = "Methylation difference",
    y = "Semantically clustered GO term",
    fill = "Genomic region"
  ) +
  facet_wrap(
    . ~ overrepresentation_status + GO_label,
    scales = "free",
    space = "fixed"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 8)),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 7.5, colour = "black", lineheight = 0.85),
    
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.text.y = element_blank(),
    
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.height = unit(0.7, "cm"),
    
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(1.2, "lines"),
    plot.margin = margin(8, 12, 8, 8)
  )

plot(ORA_exp1_methdiff_boxplot_semantic_all_terms)

ggsave(
  "ORA_exp1_methdiff_boxplot_semantic_all_terms_columns_combined_rows.png",
  ORA_exp1_methdiff_boxplot_semantic_all_terms,
  width = 42,
  height = 40,
  units = "cm",
  dpi = 600
)

## bbplot----
# -------------------------------
# Install/load bbplot
# -------------------------------

if (!requireNamespace("bbplot", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("bbc/bbplot")
}

library(bbplot)

# -------------------------------
# Professional BBC-style boxplot + points
# -------------------------------

ORA_exp1_methdiff_dotplot_semantic_all_terms <- ggplot(
  ORA_methdiff_boxplot_df,
  aes(
    x = methDiff,
    y = semantic_term_key,
    colour = region_label,
    group = region_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey45"
  ) +
  geom_point(
    position = position_jitterdodge(
      dodge.width = 0.85,
      jitter.width = 0,
      jitter.height = 0.12,
      seed = 101
    ),
    size = 1.8,
    alpha = 0.65
  ) +
  scale_colour_manual(
    values = region_cols,
    na.translate = FALSE
  ) +
  scale_y_discrete(
    labels = semantic_term_labels,
    drop = TRUE
  ) +
  labs(
    title = "Methylation differences across semantically clustered GO terms",
    subtitle = "Points show individual methylation-difference values by genomic region",
    x = "Methylation difference",
    y = NULL,
    colour = "Genomic region",
    caption = "GO terms grouped by semantic cluster description"
  ) +
  facet_wrap(
    ~ overrepresentation_status + GO_label,
    nrow = 1,
    scales = "free_y"
  ) +
  bbplot::bbc_style() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      colour = "#222222",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 13,
      colour = "#555555",
      margin = margin(b = 14)
    ),
    plot.caption = element_text(
      size = 9,
      colour = "#666666",
      hjust = 0,
      margin = margin(t = 10)
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      colour = "#222222",
      margin = margin(t = 10)
    ),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      size = 10,
      colour = "#222222"
    ),
    axis.text.y = element_text(
      size = 7.5,
      colour = "#222222",
      lineheight = 0.85
    ),
    axis.line.x = element_line(
      colour = "#222222",
      linewidth = 0.4
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    strip.background = element_rect(
      fill = "#F2F2F2",
      colour = NA
    ),
    strip.text = element_text(
      size = 10.5,
      face = "bold",
      colour = "#222222",
      margin = margin(t = 6, r = 6, b = 6, l = 6)
    ),
    
    legend.position = "bottom",
    legend.title = element_text(
      size = 10.5,
      face = "bold",
      colour = "#222222"
    ),
    legend.text = element_text(
      size = 9.5,
      colour = "#222222"
    ),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.8, "cm"),
    
    panel.grid.major.y = element_line(
      colour = "#E6E6E6",
      linewidth = 0.25
    ),
    panel.grid.major.x = element_line(
      colour = "#D9D9D9",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    panel.spacing.x = unit(1.4, "lines"),
    panel.spacing.y = unit(1.0, "lines"),
    plot.margin = margin(14, 16, 12, 14)
  ) +
  guides(
    colour = guide_legend(
      title.position = "top",
      nrow = 1,
      byrow = TRUE
    )
  )

plot(ORA_exp1_methdiff_dotplot_semantic_all_terms)

ggsave(
  "ORA_exp1_methdiff_dotplot_semantic_all_terms_bbplot_style.png",
  ORA_exp1_methdiff_dotplot_semantic_all_terms,
  width = 46,
  height = 42,
  units = "cm",
  dpi = 600,
  bg = "white"
)
plot(ORA_exp1_methdiff_boxplot_semantic_all_terms)

ggsave(
  "ORA_exp1_methdiff_boxplot_semantic_all_terms_bbplot_style.png",
  ORA_exp1_methdiff_boxplot_semantic_all_terms,
  width = 46,
  height = 34,
  units = "cm",
  dpi = 600,
  bg = "white"
)


# network cluster ----
# -------------------------------
# Semantic GO network cluster plot
# -------------------------------

# Install required packages if needed
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

if (!requireNamespace("tidygraph", quietly = TRUE)) {
  install.packages("tidygraph")
}

if (!requireNamespace("ggraph", quietly = TRUE)) {
  install.packages("ggraph")
}

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)

# -------------------------------
# 1. Summarise methylation profiles by semantic GO group
# -------------------------------
# Each semantic GO group gets a methylation profile across genomic regions.

semantic_region_profile <- ORA_methdiff_boxplot_df %>%
  group_by(
    Description_clustered,
    region_label
  ) %>%
  summarise(
    median_methDiff = median(methDiff, na.rm = TRUE),
    n_methDiff_values = n(),
    min_padj = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------
# 2. Create wide matrix for similarity calculation
# -------------------------------
# Rows = semantic GO descriptions
# Columns = genomic regions
# Values = median methylation difference
#
# Missing region values are filled with 0 so that absence from a region
# contributes neutrally to the profile.

semantic_profile_wide <- semantic_region_profile %>%
  select(
    Description_clustered,
    region_label,
    median_methDiff
  ) %>%
  pivot_wider(
    names_from = region_label,
    values_from = median_methDiff,
    values_fill = 0
  )

semantic_profile_matrix <- semantic_profile_wide %>%
  select(-Description_clustered) %>%
  as.matrix()

rownames(semantic_profile_matrix) <- semantic_profile_wide$Description_clustered

# -------------------------------
# 3. Calculate pairwise similarity between semantic GO groups
# -------------------------------
# Pearson correlation is used here to connect groups with similar
# regional methylation profiles.

semantic_similarity_matrix <- cor(
  t(semantic_profile_matrix),
  method = "pearson",
  use = "pairwise.complete.obs"
)

semantic_similarity_df <- as.data.frame(as.table(semantic_similarity_matrix)) %>%
  rename(
    from = Var1,
    to = Var2,
    similarity = Freq
  ) %>%
  mutate(
    from = as.character(from),
    to = as.character(to)
  ) %>%
  filter(
    from != to,
    !is.na(similarity)
  )

# Keep only one copy of each undirected pair
semantic_similarity_df <- semantic_similarity_df %>%
  rowwise() %>%
  mutate(
    pair_id = paste(sort(c(from, to)), collapse = "___")
  ) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)

# -------------------------------
# 4. Filter edges
# -------------------------------
# Increase similarity_cutoff for a cleaner, stricter network.
# Decrease it if too few edges are drawn.

similarity_cutoff <- 0.85

semantic_edges <- semantic_similarity_df %>%
  filter(similarity >= similarity_cutoff) %>%
  transmute(
    from,
    to,
    weight = similarity
  )

# -------------------------------
# 5. Create node metadata
# -------------------------------

semantic_nodes <- ORA_methdiff_boxplot_df %>%
  group_by(Description_clustered) %>%
  summarise(
    mean_methDiff = mean(methDiff, na.rm = TRUE),
    median_methDiff = median(methDiff, na.rm = TRUE),
    max_abs_methDiff = max(abs(methDiff), na.rm = TRUE),
    total_methDiff_values = n(),
    n_regions = n_distinct(region_label),
    min_padj = min(p.adjust, na.rm = TRUE),
    
    overrepresentation_status = if_else(
      any(overrepresentation_status == "Significantly overrepresented"),
      "Significantly overrepresented",
      "Not significantly overrepresented"
    ),
    
    GO_label = paste(sort(unique(GO_label)), collapse = " / "),
    
    .groups = "drop"
  ) %>%
  mutate(
    name = Description_clustered,
    label = str_wrap(Description_clustered, width = 28),
    overrepresentation_status = factor(
      overrepresentation_status,
      levels = c(
        "Significantly overrepresented",
        "Not significantly overrepresented"
      )
    )
  )

# Keep only nodes present in at least one retained edge
semantic_nodes_network <- semantic_nodes %>%
  filter(
    name %in% unique(c(semantic_edges$from, semantic_edges$to))
  )

semantic_edges_network <- semantic_edges %>%
  filter(
    from %in% semantic_nodes_network$name,
    to %in% semantic_nodes_network$name
  )

# -------------------------------
# 6. Build graph and detect communities
# -------------------------------

semantic_graph <- tbl_graph(
  nodes = semantic_nodes_network,
  edges = semantic_edges_network,
  directed = FALSE
) %>%
  activate(nodes) %>%
  mutate(
    community = as.factor(group_louvain(weights = weight)),
    degree = centrality_degree(),
    betweenness = centrality_betweenness()
  )

# -------------------------------
# 7. Select labels
# -------------------------------
# Label the most connected and/or significant semantic groups.

label_nodes <- semantic_graph %>%
  activate(nodes) %>%
  as_tibble() %>%
  arrange(
    min_padj,
    desc(degree),
    desc(total_methDiff_values)
  ) %>%
  slice_head(n = 25) %>%
  pull(name)

# -------------------------------
# 8. Network cluster plot
# -------------------------------

ORA_exp1_semantic_GO_network_cluster <- ggraph(
  semantic_graph,
  layout = "fr"
) +
  geom_edge_link(
    aes(
      width = weight,
      alpha = weight
    ),
    colour = "grey65",
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(
      size = total_methDiff_values,
      colour = community,
      shape = overrepresentation_status
    ),
    alpha = 0.9
  ) +
  geom_node_text(
    aes(
      label = if_else(name %in% label_nodes, label, "")
    ),
    repel = TRUE,
    size = 3,
    lineheight = 0.85,
    max.overlaps = Inf
  ) +
  scale_edge_width(
    range = c(0.2, 1.8)
  ) +
  scale_edge_alpha(
    range = c(0.15, 0.75)
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "Number of methylation\nvalues"
  ) +
  scale_shape_manual(
    values = c(
      "Significantly overrepresented" = 16,
      "Not significantly overrepresented" = 1
    ),
    name = "ORA status"
  ) +
  labs(
    title = "Network clustering of semantic GO groups",
    subtitle = paste0(
      "Edges connect semantic GO groups with similar regional methylation profiles",
      " | similarity cutoff = ",
      similarity_cutoff
    ),
    colour = "Network cluster"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      colour = "#222222",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 12,
      colour = "#555555",
      margin = margin(b = 12)
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 10.5,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9.5
    ),
    plot.margin = margin(12, 16, 12, 16)
  )

plot(ORA_exp1_semantic_GO_network_cluster)

ggsave(
  "ORA_exp1_semantic_GO_network_cluster.png",
  ORA_exp1_semantic_GO_network_cluster,
  width = 32,
  height = 24,
  units = "cm",
  dpi = 600,
  bg = "white"
)

# GOCircle plot ----

# -------------------------------
# GOCircle plot using GOplot
# Semantic GO groups as terms
# methDiff used as logFC analogue
# -------------------------------

if (!requireNamespace("GOplot", quietly = TRUE)) {
  install.packages("GOplot")
}

library(dplyr)
library(stringr)
library(tidyr)
library(GOplot)
library(ggplot2)

# -------------------------------
# 1. Prepare data for GOplot
# -------------------------------
# GOplot expects:
# terms: category, ID, term, adj_pval, genes
# genes: ID, logFC
#
# Here:
# term = semantic GO description
# genes = pseudo feature IDs representing individual methylation values
# logFC = methDiff

GOplot_input_df <- ORA_methdiff_boxplot_df %>%
  filter(
    !is.na(methDiff),
    !is.na(Description_clustered),
    !is.na(Description),
    !is.na(GO_label),
    !is.na(region_label),
    !is.na(p.adjust),
    !is.na(overrepresentation_status)
  ) %>%
  mutate(
    # GOplot usually expects BP / MF / CC-style categories
    category = recode(
      GO_label,
      "Biological process" = "BP",
      "Molecular function" = "MF",
      "Cellular component" = "CC",
      .default = GO_label
    ),
    
    # Make unique pseudo-feature IDs because individual methylation values
    # are being used as the quantitative features.
    pseudo_feature_id = paste0(
      "DMR_",
      row_number()
    ),
    
    logFC = methDiff
  )

# -------------------------------
# 2. Create gene-like quantitative table
# -------------------------------

GOplot_genes <- GOplot_input_df %>%
  distinct(
    ID = pseudo_feature_id,
    logFC
  )

# -------------------------------
# 3. Create semantic-term enrichment table
# -------------------------------
# Each semantic cluster becomes one GOCircle term.
# Original GO descriptions are collapsed underneath the semantic group.

GOplot_terms_semantic <- GOplot_input_df %>%
  group_by(
    overrepresentation_status,
    category,
    GO_label,
    Description_clustered
  ) %>%
  summarise(
    adj_pval = min(p.adjust, na.rm = TRUE),
    n_methDiff_values = n(),
    n_original_GO_terms = n_distinct(Description),
    median_methDiff = median(methDiff, na.rm = TRUE),
    
    genes = paste(unique(pseudo_feature_id), collapse = ", "),
    
    original_GO_descriptions = paste(
      sort(unique(Description)),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    overrepresentation_status,
    category,
    adj_pval,
    desc(n_methDiff_values)
  ) %>%
  mutate(
    semantic_id = paste0(
      category,
      "_S",
      stringr::str_pad(row_number(), width = 3, pad = "0")
    ),
    
    # GOplot term label.
    # Keep this readable because it appears in the GOCircle legend/table.
    term = str_trunc(
      Description_clustered,
      width = 70
    )
  ) %>%
  transmute(
    overrepresentation_status,
    GO_label,
    category,
    ID = semantic_id,
    term,
    adj_pval,
    genes,
    n_methDiff_values,
    n_original_GO_terms,
    median_methDiff,
    original_GO_descriptions
  )

# -------------------------------
# 4. Function to make one GOCircle plot
# -------------------------------
# GOCircle becomes crowded quickly, so this defaults to the top 15
# semantic groups per plot by adjusted p-value.

make_gocircle_plot <- function(
    terms_df,
    genes_df,
    status_filter,
    ontology_filter = NULL,
    top_n_terms = 15,
    plot_title = NULL,
    output_file = NULL
) {
  
  plot_terms <- terms_df %>%
    filter(
      overrepresentation_status == status_filter
    )
  
  if (!is.null(ontology_filter)) {
    plot_terms <- plot_terms %>%
      filter(
        GO_label == ontology_filter
      )
  }
  
  plot_terms <- plot_terms %>%
    arrange(
      adj_pval,
      desc(n_methDiff_values),
      desc(n_original_GO_terms)
    ) %>%
    slice_head(n = top_n_terms) %>%
    select(
      category,
      ID,
      term,
      adj_pval,
      genes
    ) %>%
    filter(
      !is.na(adj_pval),
      is.finite(adj_pval),
      adj_pval > 0,
      !is.na(genes),
      genes != ""
    )
  
  if (nrow(plot_terms) == 0) {
    stop("No valid terms available for this status / ontology combination.")
  }
  
  selected_feature_ids <- plot_terms %>%
    pull(genes) %>%
    paste(collapse = ", ") %>%
    str_split(",\\s*") %>%
    unlist() %>%
    unique()
  
  plot_genes <- genes_df %>%
    filter(
      ID %in% selected_feature_ids,
      !is.na(logFC),
      is.finite(logFC)
    )
  
  if (nrow(plot_genes) == 0) {
    stop("No matching methylation/logFC values found for the selected terms.")
  }
  
  circ <- GOplot::circle_dat(
    terms = plot_terms,
    genes = plot_genes
  )
  
  circ <- circ %>%
    filter(
      !is.na(logFC),
      is.finite(logFC),
      !is.na(adj_pval),
      is.finite(adj_pval),
      adj_pval > 0
    )
  
  if (nrow(circ) == 0) {
    stop("circle_dat() produced no valid rows after filtering.")
  }
  
  if (length(unique(circ$term)) == 0) {
    stop("No terms survived into the circle_dat() object.")
  }
  
  if (is.null(plot_title)) {
    plot_title <- paste(
      status_filter,
      ifelse(is.null(ontology_filter), "", ontology_filter)
    )
  }
  
  p <- GOplot::GOCircle(
    circ,
    title = plot_title,
    
    # Important fix:
    # use a numeric number of terms, not the synthetic semantic IDs.
    nsub = min(top_n_terms, length(unique(circ$term))),
    
    table.legend = TRUE,
    label.size = 3.5,
    label.fontface = "bold",
    lfc.col = c("#2166AC", "#B2182B"),
    zsc.col = c("#B2182B", "white", "#2166AC")
  )
  
  if (!is.null(output_file)) {
    ggsave(
      output_file,
      p,
      width = 30,
      height = 24,
      units = "cm",
      dpi = 600,
      bg = "white"
    )
  }
  
  return(p)
}
# -------------------------------
# 5. Significant semantic GO groups
# -------------------------------

GOCircle_exp1_semantic_sig <- make_gocircle_plot(
  terms_df = GOplot_terms_semantic,
  genes_df = GOplot_genes,
  status_filter = "Significantly overrepresented",
  ontology_filter = NULL,
  top_n_terms = 30,
  plot_title = "Significantly overrepresented semantic GO groups",
  output_file = "GOCircle_exp1_semantic_significant.png"
)

plot(GOCircle_exp1_semantic_sig)

# -------------------------------
# 6. Non-significant semantic GO groups
# -------------------------------

GOCircle_exp1_semantic_nonsig <- make_gocircle_plot(
  terms_df = GOplot_terms_semantic,
  genes_df = GOplot_genes,
  status_filter = "Not significantly overrepresented",
  ontology_filter = NULL,
  top_n_terms = 15,
  plot_title = "Non-significant semantic GO groups",
  output_file = "GOCircle_exp1_semantic_nonsignificant.png"
)

plot(GOCircle_exp1_semantic_nonsig)


# screw it, just a main figure table of the results to use ----












































