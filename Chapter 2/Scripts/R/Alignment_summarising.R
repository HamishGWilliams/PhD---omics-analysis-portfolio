library(readr)
library(dplyr)
library(ggplot2)
library(scales)

df <- read_tsv("C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/sequencing_and_alignment_data.txt", show_col_types = FALSE)

str(df)

plot_df <- df %>%
  mutate(
    alignment_percent = `Alignment (%)` * 100,
    x_group = factor(
      paste(Experiment, Treatment, sep = "_"),
      levels = c(
        "Acute_Control",
        "Acute_Acute",
        "Primed_Naïve",
        "Primed_Primed"
      ),
      labels = c(
        "Control",
        "Acute",
        "Naïve",
        "Primed"
      )
    ),
    Experiment = factor(Experiment, levels = c("Acute", "Primed"))
  )

str(plot_df)

group_cols <- c(
  "Control" = "#0047AB",
  "Acute"   = "#D91C3F",
  "Naïve"   = "#FFB300",
  "Primed"  = "#009E73"
)

alignment_plot <- ggplot(plot_df, aes(x = x_group, y = alignment_percent, fill = x_group)) +
  geom_violin(
    width = 0.75,
    alpha = 0.60,
    colour = "grey30"
  ) +
  geom_jitter(
    width = 0.12,
    size = 3,
    alpha = 0.75,
    colour = "black"
  ) +
  geom_vline(xintercept = 2.5, linewidth = 0.4, colour = "black") +
  annotate("text", x = 1.5, y = max(plot_df$alignment_percent) + 1, label = "Acute") +
  annotate("text", x = 3.5, y = max(plot_df$alignment_percent) + 1, label = "Primed") +
  labs(
    x = NULL,
    y = "Alignment (%)"
  ) +
  scale_fill_manual(values = group_cols) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.margin = margin(10, 10, 20, 10)
  )

alignment_plot

ggsave(
  filename = file.path("C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Results/alignment_violin.png"),
  plot = alignment_plot,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)
           

# Wilcoxons paired test on alignment efficiencies: ----
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

df <- read_tsv(
  "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/Chapter 2/Data/sequencing_and_alignment_data.txt",
  show_col_types = FALSE
)

plot_df <- df %>%
  mutate(
    alignment_percent = `Alignment (%)` * 100,
    x_group = factor(
      paste(Experiment, Treatment, sep = "_"),
      levels = c(
        "Acute_Control",
        "Acute_Acute",
        "Primed_Naïve",
        "Primed_Primed"
      ),
      labels = c(
        "Control",
        "Acute",
        "Naïve",
        "Primed"
      )
    ),
    Experiment = factor(Experiment, levels = c("Acute", "Primed"))
  )

pair_map_exp1 <- c(
  "3"  = "5",
  "9"  = "15",
  "13" = "9",
  "15" = "15",
  "16" = "2",
  "17" = "5",
  "22" = "9",
  "27" = "2"
)

pair_map_exp2 <- c(
  "6"  = "19",
  "7"  = "8",
  "8"  = "20",
  "10" = "6",
  "14" = "8",
  "20" = "19",
  "21" = "12",
  "11" = "12",
  "18" = "20",
  "2"  = "6"
)

paired_df <- plot_df %>%
  mutate(
    sample_id = as.character(sample_id),
    pair_id = case_when(
      Experiment == "Acute"  ~ unname(pair_map_exp1[sample_id]),
      Experiment == "Primed" ~ unname(pair_map_exp2[sample_id]),
      TRUE ~ NA_character_
    ),
    comparison_group = case_when(
      Experiment == "Acute"  & Treatment == "Control" ~ "Control",
      Experiment == "Acute"  & Treatment == "Acute"   ~ "Acute",
      Experiment == "Primed" & Treatment == "Naïve"   ~ "Naïve",
      Experiment == "Primed" & Treatment == "Primed"  ~ "Primed",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pair_id), !is.na(comparison_group))

paired_df

pair_check <- paired_df %>%
  count(Experiment, pair_id, comparison_group) %>%
  tidyr::pivot_wider(
    names_from = comparison_group,
    values_from = n,
    values_fill = 0
  )

pair_check

paired_wide <- paired_df %>%
  select(Experiment, pair_id, comparison_group, alignment_percent) %>%
  pivot_wider(
    names_from = comparison_group,
    values_from = alignment_percent
  )

paired_wide

wilcox_results <- bind_rows(
  paired_wide %>%
    filter(Experiment == "Acute") %>%
    filter(!is.na(Control), !is.na(Acute)) %>%
    summarise(
      Experiment = "Acute",
      n_pairs = n(),
      statistic = unname(wilcox.test(Control, Acute, paired = TRUE, exact = FALSE)$statistic),
      p_value = wilcox.test(Control, Acute, paired = TRUE, exact = FALSE)$p.value
    ),
  
  paired_wide %>%
    filter(Experiment == "Primed") %>%
    filter(!is.na(Naïve), !is.na(Primed)) %>%
    summarise(
      Experiment = "Primed",
      n_pairs = n(),
      statistic = unname(wilcox.test(Naïve, Primed, paired = TRUE, exact = FALSE)$statistic),
      p_value = wilcox.test(Naïve, Primed, paired = TRUE, exact = FALSE)$p.value
    )
) %>%
  mutate(
    method = "Wilcoxon signed-rank test",
    p_adj_BH = p.adjust(p_value, method = "BH"),
    p_signif = case_when(
      p_value <= 0.0001 ~ "****",
      p_value <= 0.001  ~ "***",
      p_value <= 0.01   ~ "**",
      p_value <= 0.05   ~ "*",
      TRUE ~ "ns"
    )
  )

wilcox_results

