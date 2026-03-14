# Single R script for making the upset plot comparing all contrasts and methods together

# Load packages as usual ----
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

# load in merged_df with all the data needeed (nice)
data <- read.csv("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Results/Multi_Stressor_all_Differential_Expression_Analysis_results.csv",
                 header = T,
                 sep = ",")

str()

# thresholds
alpha   <- 0.05      # FDR threshold
lfc_thr <- NA_real_  # set e.g. 0.58 for ~1.5x threshold on WALD only; keep NA to ignore

# filter data in results to only DEGs
sets <- data %>%
  transmute(
    gene,
    # LRT values
    `LRT: D`            = !is.na(padj_diesel_only_LRT)   & padj_diesel_only_LRT   < alpha,
    `LRT: S`          = !is.na(padj_salinity_only_LRT) & padj_salinity_only_LRT < alpha,
    `LRT: D + S + D*S`  = !is.na(padj.x)      & padj.x      < alpha,
    `LRT: D*S`  = !is.na(padj.y)      & padj.y      < alpha,
    
    # Wald values
    `Wald: D`           = !is.na(padj_diesel_only_wald)   & padj_diesel_only_wald   < alpha,
    `Wald: S`         = !is.na(padj_salinity_only_wald) & padj_salinity_only_wald < alpha,
    `Wald: +D`           = !is.na(padj_diesel_added_wald)   & padj_diesel_added_wald   < alpha,
    `Wald: +S`         = !is.na(padj_salinity_added_wald) & padj_salinity_added_wald < alpha,
    `Wald: D + S + D*S` = !is.na(padj_combined_wald)      & padj_combined_wald      < alpha,
    `Wald: D*S` = !is.na(padj_interactive_only_wald)      & padj_interactive_only_wald      < alpha
  )

sets_wald <- data %>%
  transmute(
    gene,
    # Wald values
    `Diesel`           = !is.na(padj_diesel_only_wald)   & padj_diesel_only_wald   < alpha,
    `Salinity`         = !is.na(padj_salinity_only_wald) & padj_salinity_only_wald < alpha,
    `+ Diesel`           = !is.na(padj_diesel_added_wald)   & padj_diesel_added_wald   < alpha,
    `+ Salinity`         = !is.na(padj_salinity_added_wald) & padj_salinity_added_wald < alpha,
    `D + S + D*S` = !is.na(padj_combined_wald)      & padj_combined_wald      < alpha,
    `D*S` = !is.na(padj_interactive_only_wald)      & padj_interactive_only_wald      < alpha
  )


# Make upset plot ----
m1 <- make_comb_mat(sets_wald)

# remove the empty (degree 0) intersection; also drop zero-size combos just in case
m1 <- m1[ comb_degree(m1) > 0 & comb_size(m1) > 0 ]

png("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Figures/Main Figures/upset_plot.png", width = 4000, height = 2000, res = 300)

UpSet(m1,
      set_order = c("Diesel",
                    "Salinity",
                    "+ Diesel",
                    "+ Salinity",
                    "D + S + D*S",
                    "D*S"
                    ),
      
      comb_order = order(comb_size(m1), decreasing = TRUE),
      
      top_annotation = upset_top_annotation(
        m1,
        add_numbers   = TRUE,
        numbers_gp    = gpar(cex = 0.7, fontface = "bold"),
        numbers_offset= unit(1, "mm"),
        numbers_rot   = 0,
        gp            = gpar(fill = "grey50", col = "grey30"),
        bar_width     = unit(6, "mm"),
        ylim          = c(0, max(comb_size(m1))),
        axis_param    = list(
          at      = pretty(c(0, max(comb_size(m1)))),
          labels  = pretty(c(0, max(comb_size(m1)))),
          gp      = gpar(cex = 0.75)
        ),
        height        = unit(40, "mm")
      ),
      
      right_annotation = upset_right_annotation(
        m1,
        add_numbers = TRUE,
        numbers_gp  = gpar(cex = 0.7),
        gp          = gpar(fill = "grey75", col = "grey30"),
        width       = unit(4, "cm"),  # 1) Make the set-size barplot region wider (longer bars)
        bar_width   = unit(6, "mm")
      ),
      
      comb_col = "black", 
      #bg_col = rep(c("white", "grey97"), 3),
      bg_col = c("#cc4551","#cc9593", "#f7e9a1", "#98c2d9", "#f7e9a1","#98c2d9"), 
      bg_pt_col = "white",
      pt_size = unit(3, "mm"),
      lwd = 1.2,
      row_names_gp = gpar(fontsize = 12)
      )


dev.off()


# Summary table of DEG directions
summary_table <- data %>%
  rename(
    padj_full_model_LRT = padj.x,
    padj_interactive_only_LRT = padj.y
  ) %>%
  pivot_longer(
    cols = -gene,
    names_to = c("metric", "comparison"),
    names_pattern = "^(LFC|padj)_(.*)$",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  ) %>%
  mutate(
    direction = case_when(
      padj < 0.05 & LFC > 0 ~ "Upregulated",
      padj < 0.05 & LFC < 0 ~ "Downregulated",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(direction)) %>%
  count(comparison, direction, name = "n_genes") %>%
  pivot_wider(
    names_from = direction,
    values_from = n_genes,
    values_fill = 0
  ) %>%
  mutate(total_DE = Upregulated + Downregulated)

summary_table

# comparison                  Downregulated   Upregulated   total_DE
# 
# 1 combined_wald             4250            3360          7610
# 2 diesel_added_wald         999             626           1625
# 3 diesel_only_wald          3               4             7
# 4 full_model_LRT            5301            5543          10844
# 5 interactive_only_LRT      37              20            57
# 6 interactive_only_wald     70              55            125
# 7 salinity_added_wald       4067            3576          7643
# 8 salinity_only_LRT         2999            2440          5439
# 9 salinity_only_wald        3540            2872          6412