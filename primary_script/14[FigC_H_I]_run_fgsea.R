# ----Libraries----
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(GSEABase)

library(tidyverse)
library(pheatmap)

library(ggvenn)
library(cowplot)

library(tidytext)   

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("lag", "dplyr")

library(here)


# ----Load data----
# DGE SD1-3 (cMono, nMono, MPA): written by Script 9
mast_all_monos_sd1_sd3_deg <- read_rds(here::here("intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds"))

# DGE MPAs vs cMonos without platelet genes: written by Script 9_13
mast_mpa_cmono_sd1_wo_platelet_genes_deg <- read_csv(here::here("intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.csv"))
mast_mpa_cmono_sd3_wo_platelet_genes_deg <- read_csv(here::here("intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.csv"))

# DGE MPA real vs sim without platelet genes: written by Script 11
mast_sim_dblts_without_platelet_genes_deg <- read_rds(here::here("intermediate/mast/mpa_real_sim/pbmc_mpa_sim_real_no_platelet_genes_deg_6FEB2025.rds"))

# ----Rank for fgsea----

mpa_sd1_sd3_rank <- mast_all_monos_sd1_sd3_deg %>%
  filter(cluster == "M-platelet") %>%          # Filter for M-platelet cluster
  arrange(desc(logFC)) %>%                    # Arrange by logFC (descending)
  select(primerid, logFC) %>%                 # Keep only needed columns
  deframe()                                   # Convert to named numeric vector (primerid = name)

cmono_sd1_sd3_rank <- mast_all_monos_sd1_sd3_deg %>%
  filter(cluster == "CD14 Mono") %>%          # Filter for M-platelet cluster
  arrange(desc(logFC)) %>%                    # Arrange by logFC (descending)
  select(primerid, logFC) %>%                 # Keep only needed columns
  deframe()                                   # Convert to named numeric vector (primerid = name)

nmono_sd1_sd3_rank <- mast_all_monos_sd1_sd3_deg %>%
  filter(cluster == "CD16 Mono") %>%          # Filter for M-platelet cluster
  arrange(desc(logFC)) %>%                    # Arrange by logFC (descending)
  select(primerid, logFC) %>%                 # Keep only needed columns
  deframe()                                   # Convert to named numeric vector (primerid = name)

mpa_cmono_sd1_wo_plt_rank <- mast_mpa_cmono_sd1_wo_platelet_genes_deg %>%
  arrange(desc(logFC)) %>%                    # Arrange by logFC (descending)
  select(primerid, logFC) %>%                 # Keep only needed columns
  deframe()                                   # Convert to named numeric vector (primerid = name)

mpa_cmono_sd3_wo_plt_rank <- mast_mpa_cmono_sd3_wo_platelet_genes_deg %>%
  arrange(desc(logFC)) %>%                    # Arrange by logFC (descending)
  select(primerid, logFC) %>%                 # Keep only needed columns
  deframe()                                   # Convert to named numeric vector (primerid = name)

sim_dble_wo_plt_rank <- mast_sim_dblts_without_platelet_genes_deg %>%
  arrange(desc(logFC)) %>%                    # Arrange by logFC (descending)
  select(primerid, logFC) %>%                 # Keep only needed columns
  deframe()                                   # Convert to named numeric vector (primerid = name)


# ----Pick gene set----

# Hallmark
gene_set <- GSEABase::getGmt(here::here("primary_dependents/gene_set/H_Hallmark_modules_gsea.gmt"))

# where .gmt file can be any gmt file; for other msigdb gmt files, see: https://github.com/jsim91/seutools/tree/master/inst/extdata
gene_ids <- GSEABase::geneIds(gene_set)

# ----Set seed----
set.seed(123)

# ----Run fgsea----
fgseaRes_mpa_sd1_sd3 <- fgsea(pathways = gene_ids, 
                  stats    = mpa_sd1_sd3_rank,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

fgseaRes_cmono_sd1_sd3 <- fgsea(pathways = gene_ids, 
                              stats    = cmono_sd1_sd3_rank,
                              eps      = 0.0,
                              minSize  = 15,
                              maxSize  = 500)

fgseaRes_nmono_sd1_sd3 <- fgsea(pathways = gene_ids, 
                                stats    = nmono_sd1_sd3_rank,
                                eps      = 0.0,
                                minSize  = 15,
                                maxSize  = 500)


fgseaRes_mpa_cmono_sd1_wo_plt <- fgsea(pathways = gene_ids, 
                              stats    = mpa_cmono_sd1_wo_plt_rank,
                              eps      = 0.0,
                              minSize  = 15,
                              maxSize  = 500)

fgseaRes_mpa_cmono_sd3_wo_plt <- fgsea(pathways = gene_ids, 
                                stats    = mpa_cmono_sd3_wo_plt_rank,
                                eps      = 0.0,
                                minSize  = 15,
                                maxSize  = 500)

fgseaRes_sim_dble_wo_plt <- fgsea(pathways = gene_ids, 
                                 stats    = sim_dble_wo_plt_rank,
                                 eps      = 0.0,
                                 minSize  = 15,
                                 maxSize  = 500)

# ----Sim Dblt Figure----
# ---- FGSEA Bar Plot Function (NES-based) 
visualize_fgsea_barplots <- function(fgsea_results, 
                                     output_dir = "fgsea_barplots_nes", 
                                     significance_cutoff = 0.05) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Filter significant pathways
  significant_results <- fgsea_results %>%
    filter(padj < significance_cutoff) %>%
    arrange(NES)  # Sort by NES (low to high or reverse later if preferred)
  
  # Check if there are significant results to plot
  if (nrow(significant_results) == 0) {
    message("No significant pathways found to plot.")
    return(NULL)
  }
  
  # Prepare data for plotting
  gsea_df <- significant_results %>%
    mutate(pathway = factor(pathway, levels = pathway),  # Maintain order
           log_padj = -log10(padj))                      # Transform padj for color
  
  # Plot NES with color scaled by -log10(padj)
  bar_plot <- ggplot(gsea_df, aes(x = NES, y = pathway, fill = log_padj)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#a6cee3", high = "#1f78b4", name = "-log10(adj p)") +
    labs(title = "FGSEA - NES Bar Plot (Significant Pathways)",
         x = "Normalized Enrichment Score (NES)",
         y = "Pathway") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6),  # Adjust as needed
          axis.text.x = element_text(size = 8),
          plot.title = element_text(size = 10, face = "bold"))
  
  # Save the plot
  ggsave(filename = file.path(output_dir, "fgsea_NES_barplot.png"),
         plot = bar_plot,
         width = 10,
         height = max(4, nrow(gsea_df) * 0.3), #Dynamically adjust height
         limitsize = FALSE,
         dpi = 300)
  
  return(bar_plot)
}


# Run the visualization function
sim_dble_wo_plt <- visualize_fgsea_barplots(fgseaRes_sim_dble_wo_plt)

dir.create(here::here("output/figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(here::here("output/figures/sim_dblt_wo_plt_gsea.pdf"), sim_dble_wo_plt, width = 8, height = 6)

# ----Combined Monos Figure----

# sd1_sd3
fgsea_all_monos_list_named <- list(
  cMono = fgseaRes_cmono_sd1_sd3,
  nMono = fgseaRes_nmono_sd1_sd3,
  MPA   = fgseaRes_mpa_sd1_sd3
)

combined_fgsea <- fgsea_all_monos_list_named %>%
  imap_dfr(~ .x %>% 
             mutate(cell_type = .y) %>%
             select(pathway, NES, padj, cell_type))

# Optional - filter for pathways significant somewhere
combined_fgsea_filtered <- combined_fgsea %>% 
  filter(padj < 0.05)

# Optional - find common pathways across all cell types
common_pathways <- combined_fgsea_filtered %>%
  distinct(pathway, cell_type) %>%
  count(pathway) %>%
  filter(n == length(fgsea_all_monos_list_named)) %>%
  pull(pathway)

combined_fgsea_common_filtered <- combined_fgsea_filtered %>%
  filter(pathway %in% common_pathways)

pathway_order <- combined_fgsea_common_filtered %>%
  filter(cell_type == "cMono") %>%
  arrange(desc(NES)) %>%
  pull(pathway)

combined_sd1_sd3 <- combined_fgsea_common_filtered %>%
  mutate(pathway = factor(pathway, levels = rev(pathway_order))) %>%
  ggplot(aes(x = NES, y = pathway, fill = -log10(padj))) +
  geom_col() +
  geom_text(aes(label = round(NES, 2)),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5) +
  facet_wrap(~cell_type, scales = "free_x") +
  theme_minimal() +
  labs(title = "Significant Pathways (Baseline to 12-Weeks)",
       x = "Normalized Enrichment Score (NES)",
       y = NULL) + 
  scale_fill_gradient(low = "#a6cee3", high = "#1f78b4", name = "-log10(adj p)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.title.position = "plot",
    strip.text = element_text(size = 12, face = "bold"),      # Facet labels
    axis.title.x = element_text(size = 12, face = "bold"),    # X-axis title
    axis.title.y = element_text(size = 12, face = "bold")     # Y-axis title
  )


ggsave(here::here("output/figures/combined_sd1_sd3_gsea.pdf"), combined_sd1_sd3, width = 10, height = 6)


# ----MPA vs cMono Figure----

visualize_combined_fgsea_barplots <- function(fgsea_input_list, 
                                              significance_cutoff = 0.05,
                                              output_file = "fgsea_combined_NES_barplot.png",
                                              output_dir = "fgsea_barplots_combined") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # 1) Combine & filter
  combined_fgsea <- imap_dfr(fgsea_input_list, ~ .x %>%
                               filter(padj < significance_cutoff) %>%
                               mutate(timepoint = .y,
                                      log_padj  = -log10(padj)))
  
  if (nrow(combined_fgsea) == 0) {
    message("No significant pathways found across timepoints.")
    return(NULL)
  }
  
  # 2) Reorder pathway *within* each timepoint by NES
  combined_fgsea <- combined_fgsea %>%
    mutate(pathway = reorder_within(pathway, NES, timepoint))
  
  # 3) Plot with facet‐aware y‐ordering
  bar_plot <- ggplot(combined_fgsea, aes(x = NES, y = pathway, fill = log_padj)) +
    geom_col() +
    facet_wrap(~timepoint, scales = "free_y") +
    scale_y_reordered() +
    scale_fill_gradient(low = "#F6C141", high = "#D55E00", name = "-log10(adj p)") +
    labs(title = "FGSEA – NES Bar Plot (Significant Pathways)",
         x     = "Normalized Enrichment Score (NES)",
         y     = "Pathway") +
    theme_minimal() +
    theme(
      axis.text.y  = element_text(size = 6),
      axis.text.x  = element_text(size = 8),
      strip.text   = element_text(size = 10, face = "bold"),
      plot.title   = element_text(size = 12, face = "bold")
    )
  
  # 4) Save
  ggsave(filename  = file.path(output_dir, output_file),
         plot      = bar_plot,
         width     = 12,
         height    = max(5, nrow(combined_fgsea) * 0.2),
         limitsize = FALSE,
         dpi       = 300)
  
  return(bar_plot)
}


mpa_cmono <- visualize_combined_fgsea_barplots(
  fgsea_input_list = list(
    SD1 = fgseaRes_mpa_cmono_sd1_wo_plt,
    SD3 = fgseaRes_mpa_cmono_sd3_wo_plt
  )
)
ggsave(here::here("output/figures/mpa_cmono_gsea.pdf"),
       plot   = mpa_cmono,
       width  = 12,
       height = 6)
