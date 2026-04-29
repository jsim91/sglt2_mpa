# ============================================================================
# SADIE 10x MM-148 Cell Count Analysis and Visualization
# Author: M Mashayekhi
# Purpose: Analyze and visualize cluster frequencies across study timepoints
# ============================================================================

# Analyze and visualize cluster frequencies across study timepoints
# Load libraries ----
library(Hmisc)        # also loads ggplot2, lattice, survival, Formula
library(tidyverse)
library(knitr)
library(tables)       # Table generation
library(htmlTable)
library(ggiraph)
library(vtree)
library(ggpubr)
library(sfsmisc)
library(trelliscopejs)
library(stringr)

# Load data ----
# Load myeloid cell counts by type and sample (output from script 12)
mm148_myeloid_counts <- read_csv(here::here("intermediate/pbmc/anndata_elements/mm_myeloid_named_cluster_cell_counts.csv"))

# Dataframe setup ----

# Read myeloid counts (first column contains subject_id_study_day identifiers from Script 12)
mm148_myeloid_raw <- read_csv(here::here("intermediate/pbmc/anndata_elements/mm_myeloid_named_cluster_cell_counts.csv"))
col_rename <- colnames(mm148_myeloid_raw)[1]

# Extract subject ID and study day from first column (format: subject_id_study_day)
df_myeloid_all <- mm148_myeloid_raw %>%
  rename(pt_day_id = !!col_rename) %>%
  separate(pt_day_id, into = c("subject_id", "study_day"), sep = "_", remove = TRUE) %>%
  select(subject_id, study_day, everything())

# Pull column names (myeloid cell types)
col_names_myeloid <- colnames(df_myeloid_all)
col_names_myeloid <- col_names_myeloid[!col_names_myeloid %in% c("subject_id", "study_day")]

# Standardize study day factor levels
df_myeloid_all <- df_myeloid_all %>%
  mutate(
    study_day = factor(study_day,
                       levels = c("Baseline", "2-Weeks", "12-Weeks")))

# Pivot to long format
df_myeloid_all_long <- df_myeloid_all %>%
  pivot_longer(cols = all_of(col_names_myeloid),
               names_to = "cluster",
               values_to = "value")

# Calculate frequencies across all myeloid cells per sample
df_myeloid_freq_long <- df_myeloid_all_long %>%
  group_by(subject_id, study_day) %>%
  mutate(total_myeloid = sum(value, na.rm = TRUE)) %>%
  mutate(frequency = (value / total_myeloid) * 100) %>%
  ungroup() %>%
  select(subject_id, study_day, cluster, frequency) %>%
  rename(value = frequency)

# Set theme and graphics ----

# Set ggplot theme
theme_set(theme_bw() +
            theme(axis.text.x = (element_text(size = 16,
                                            hjust = 1)),
                  axis.title.x = (element_text(size = 20)),
                  axis.text.y = (element_text(size = 16)),
                  axis.title.y = (element_text(size = 20)),
                  plot.caption = element_text(face = "italic",
                                              size = 6),
                  legend.position = "none",
                  legend.background = element_rect(
                    fill = "transparent"),
                  legend.text = element_text(size = 12),
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90",
                                                    linetype = 2),
                  plot.title = element_text(size = 12)
          ))

options(grType='plotly') # for certain graphics functions

study_day_colors <- c("Baseline" = "#1f78b4",
               "2-Weeks" = "#1f78b4",
               "12-Weeks" = "#1f78b4")

# Graph functions ----

plot_violin <- function(df, x_label = NULL, y_label = NULL, pd = 0.12, ...) {
  ggplot(df, aes(x = study_day, y = value)) +
    geom_violin(aes(color = study_day, fill = study_day), trim = FALSE,
                alpha  = 0.05) +
    geom_line(aes(group = study_id), size = 0.2) + # Adjusted line thickness
    geom_point(aes(fill = study_day), pch = 21, alpha = 1, stroke = 0.01, size = 7.5) +
    geom_boxplot(aes(color = study_day, fill = study_day), outlier.shape = NA, width = 0.2, fatten = 1, coef = 0, lwd = 0.5, alpha = 0.1) + # Removed whiskers
    scale_color_manual(values = study_day_colors) +
    scale_fill_manual(values = study_day_colors) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = ggtext::element_markdown(size = 28),
      axis.text.x = element_text(size = 28),
      axis.text.y = element_text(size = 28)
    ) +
    labs(
      x = x_label, y = y_label,
      shape = "", fill = "", linetype = "", color = ""
    )
}

# Generate figures ----

# Figure 1: MPA cell frequency
pl_violin_mpa <- plot_violin(
  df_myeloid_freq_long %>%
    rename(study_id = subject_id) %>%
    filter(cluster == "M-platelet"),
  variable_name = "cluster",
  y_label = "MPA<br>(% of myeloid cells)")

print(pl_violin_mpa)

# Figure 2: cMono cell frequency
pl_violin_cmono <- plot_violin(
  df_myeloid_freq_long %>%
    rename(study_id = subject_id) %>%
    filter(cluster == "CD14 Mono"),
  variable_name = "cluster",
  y_label = "cMono<br>(% of myeloid cells)")

print(pl_violin_cmono)

# Figure 3: Combined MPA and cMono
pl_mpa_cmono_legend <- cowplot::get_legend(pl_violin_mpa)

plot_mpa_cmono <-
  cowplot::plot_grid( pl_violin_mpa + theme(legend.position = "none"),
                      pl_violin_cmono + theme(legend.position = "none"),
                      nrow = 1)

cowplot::plot_grid( plot_mpa_cmono,
                    pl_mpa_cmono_legend,
                    ncol = 1, rel_heights = c(1,.1))

# Optional: Save figure
# ggsave(plot_mpa_cmono,
#        filename = here::here("figures/mm148/mpa_cmono_frequency.pdf")
# )
