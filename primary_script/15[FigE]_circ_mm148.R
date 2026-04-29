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
# Load the most recent data (output from script 12)
mm148_pbmc_all <- read_csv(here::here("intermediate/pbmc/anndata_elements/mm_pbmc_named_cluster_cell_counts.csv"))

# Dataframe setup ----

# Pull out ID from pt_day_id column
df_pbmc_all <- mm148_pbmc_all %>%
  mutate(subject_id = str_to_lower(str_remove_all(
    str_extract(pt_day_id, "(Sadie-\\d{3})"), "-")),
    study_day = str_to_lower(str_extract(pt_day_id, "SD\\d"))) %>%
  select(-pt_day_id) %>%
  select(subject_id, study_day, everything())

# Pull column names
col_names_pbmc <- colnames(df_pbmc_all)
col_names_pbmc <- col_names_pbmc[!col_names_pbmc %in% c("subject_id", "study_day" )]

# Modify study days
df_pbmc_all <- df_pbmc_all %>%
  mutate(
    study_day = case_when(
      study_day %in%
        c("sd1") ~
        factor("Baseline", levels =
                 c("Baseline", "2-Weeks",
                   "12-Weeks")),
      study_day %in%
        c("sd2") ~
        factor("2-Weeks", levels =
                 c("Baseline", "2-Weeks",
                   "12-Weeks")),
      study_day %in%
        c("sd3") ~
        factor("12-Weeks", levels =
                 c("Baseline", "2-Weeks",
                   "12-Weeks"))))

# Make long DF
df_pbmc_all_long <- df_pbmc_all %>%
  pivot_longer(cols = col_names_pbmc,
               names_to = "cluster",
               values_to = "value")

# Make it frequency of subset of cells

# Define subsets based on cell type groupings
b_cell_clusters <- c('B-platelet', 'Bint', 'Bmem', 'Bnaive', 'Plasmablast')
t_nk_cell_clusters <- c('CD16 NK', 'CD8mem', 'NKT', 'Treg', 'Prolif NK', 'T-platelet', 'dnT', 'gdT/MAIT', 'CD4mem', 'CD4naive', 'CD56-bright NK', 'CD8naive', 'B-NK')
myeloid_cell_clusters <- c('CD14 Mono', 'CD16 Mono', 'pDC', 'M-platelet', 'cDC1', 'cDC2')

# Define second set of subsets to include platelets with myeloid
myeloid_cell_w_plt_cluster <- c(myeloid_cell_clusters, 'Platelet')

# Assign a subset label to each cluster
df_pbmc_all_long <- df_pbmc_all_long %>%
  mutate(
    subset = case_when(
      cluster %in% b_cell_clusters ~ "B cells",
      cluster %in% t_nk_cell_clusters ~ "T/NK cells",
      cluster %in% myeloid_cell_clusters ~ "Myeloid cells",
      TRUE ~ "Other"  # Default case for unexpected clusters
    ),
    subset_2 = case_when(
      cluster %in% myeloid_cell_w_plt_cluster ~ "Myeloid w plt",
      TRUE ~ "Other"  # Default case for unexpected clusters
    )
  )

# Calculate subset totals and frequencies
df_pbmc_subset1_long <- df_pbmc_all_long %>%
  group_by(subject_id, study_day, subset) %>%
  mutate(subset_total = sum(value, na.rm = TRUE)) %>%
  mutate(frequency = (value / subset_total) *100) %>%
  ungroup() %>%
  select(-subset_total, -subset, -value, -subset_2) %>%   # Remove intermediate totals if not needed
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
  df_pbmc_subset1_long %>%
    rename(study_id = subject_id) %>%
    filter(cluster == "M-platelet"),
  variable_name = "cluster",
  y_label = "MPA<br>(% of myeloid cells)")

print(pl_violin_mpa)

# Figure 2: cMono cell frequency
pl_violin_cmono <- plot_violin(
  df_pbmc_subset1_long %>%
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
