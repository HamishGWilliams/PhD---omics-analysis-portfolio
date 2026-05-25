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

# filter data in results to only DEGs - Not including LRT results anymore for upset plot
sets <- data %>%
  transmute(
    gene,
    # LRT values
    `LRT: D`            = !is.na(padj_diesel_only_LRT)   & padj_diesel_only_LRT   < alpha,
    `LRT: S`            = !is.na(padj_salinity_only_LRT) & padj_salinity_only_LRT < alpha,
    `LRT: D + S + D*S`  = !is.na(padj.x)      & padj.x      < alpha,
    `LRT: D*S`          = !is.na(padj.y)      & padj.y      < alpha,
    
    # Wald values
    `Wald: D`            = !is.na(padj_diesel_only_wald)   & padj_diesel_only_wald   < alpha,
    `Wald: S`            = !is.na(padj_salinity_only_wald) & padj_salinity_only_wald < alpha,
    `Wald: +D`           = !is.na(padj_diesel_added_wald)   & padj_diesel_added_wald   < alpha,
    `Wald: +S`           = !is.na(padj_salinity_added_wald) & padj_salinity_added_wald < alpha,
    `Wald: D + S + D*S`  = !is.na(padj_combined_wald)      & padj_combined_wald      < alpha,
    `Wald: D*S`          = !is.na(padj_interactive_only_wald)      & padj_interactive_only_wald      < alpha
  )

# make upset matrix for later
m2 <- make_comb_mat(sets)

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
# Load environment image to manipulate plot whenever:
load("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/RData/upset_plot.RData")


# Start image save
png("C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/Figures/Main Figures/upset_plot.png", 
    width = 3250, 
    height = 2000, 
    res = 300)


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
        numbers_gp    = gpar(cex = 0.6), # , fontface = "bold"
        numbers_offset= unit(1, "mm"),
        numbers_rot   = 0,
        gp            = gpar(fill = "grey50", col = "grey30"),
        bar_width     = unit(5, "mm"),
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


# Summary table of DEG directions ----
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

# Intersections Table: ----
library(dplyr)
library(tibble)
library(ComplexHeatmap)

comb_codes <- comb_name(m2)
set_names  <- set_name(m2)

intersection_table <- tibble(
  comb_code = comb_codes,
  n_comparisons = comb_degree(m2),
  n_genes = comb_size(m2),
  intersection = vapply(comb_codes, function(code) {
    hits <- set_names[as.logical(as.integer(strsplit(code, "")[[1]]))]
    paste(hits, collapse = " & ")
  }, character(1))
) %>%
  filter(n_genes > 0) %>%
  select(intersection, n_comparisons, n_genes) %>%
  arrange(desc(n_genes), desc(n_comparisons))

intersection_table_with_total <- bind_rows(
  intersection_table,
  intersection_table %>%
    summarise(
      intersection = "TOTAL",
      n_comparisons = NA_integer_,
      n_genes = sum(n_genes)
    )
)

# Make custom function to find intersections between specific groups ----
library(dplyr)
library(tidyr)

get_specific_intersection <- function(sets_df,
                                      target_sets,
                                      gene_col = "gene",
                                      exact = FALSE,
                                      return_genes = TRUE) {
  
  # Check inputs
  missing_sets <- setdiff(target_sets, names(sets_df))
  if (length(missing_sets) > 0) {
    stop("These set names were not found in the data: ",
         paste(missing_sets, collapse = ", "))
  }
  
  comparison_cols <- setdiff(names(sets_df), gene_col)
  other_sets <- setdiff(comparison_cols, target_sets)
  
  # Make sure membership columns are logical and NAs become FALSE
  dat <- sets_df %>%
    mutate(
      across(all_of(comparison_cols), ~ replace_na(as.logical(.x), FALSE))
    )
  
  # Keep genes in all requested sets
  dat <- dat %>%
    filter(if_all(all_of(target_sets), ~ .x))
  
  # If exact = TRUE, remove genes that are also in any other set
  if (exact) {
    dat <- dat %>%
      filter(if_all(all_of(other_sets), ~ !.x))
  }
  
  # Add a readable label
  dat <- dat %>%
    mutate(
      requested_intersection = paste(target_sets, collapse = " & ")
    )
  
  if (return_genes) {
    out <- list(
      intersection = paste(target_sets, collapse = " & "),
      exact = exact,
      n_genes = nrow(dat),
      genes = dat[[gene_col]],
      data = dat
    )
  } else {
    out <- tibble(
      intersection = paste(target_sets, collapse = " & "),
      exact = exact,
      n_genes = nrow(dat)
    )
  }
  
  return(out)
}

# Use function
  # Names:
    # # LRT values
    # `LRT: D`           
    # `LRT: S`           
    # `LRT: D + S + D*S` 
    # `LRT: D*S`         
    # # Wald values
    # `Wald: D`          
    # `Wald: S`          
    # `Wald: +D`         
    # `Wald: +S`         
    # `Wald: D + S + D*S`
    # `Wald: D*S`        

res_inc <- get_specific_intersection(
  sets_df = sets,
  target_sets = c("LRT: D*S", "Wald: D*S"),
  exact = FALSE, # ? only include genes unique to this contrast?
  return_genes = FALSE # ? Capture gene names too?
  )

res_inc$n_genes

# multiple intersections at once: ----
summarise_requested_intersections <- function(sets_df, requests, gene_col = "gene", exact = FALSE) {
  bind_rows(lapply(requests, function(x) {
    get_specific_intersection(
      sets_df = sets_df,
      target_sets = x,
      gene_col = gene_col,
      exact = exact,
      return_genes = FALSE
    )
  }))
}

requests <- list(
  c("LRT: D", "Wald: D"),
  c("LRT: S", "Wald: S"),
  c("Wald: +D", "Wald: D + S + D*S")
)

summarise_requested_intersections(sets, requests, exact = FALSE)


# save to load later ----
save.image(file = "C:/Users/hamis/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 3/RData/upset_plot.RData")
