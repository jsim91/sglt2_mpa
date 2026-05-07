rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(matrixStats)
  library(ggdendro)
  library(patchwork)
})

resolve_repo_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "primary_script")) && dir.exists(file.path(current, "intermediate"))) {
      return(current)
    }

    candidate <- file.path(current, "technical_review", "submission_partition_draft")
    if (dir.exists(file.path(candidate, "primary_script")) && dir.exists(file.path(candidate, "intermediate"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop("Could not locate submission_partition_draft root.")
}

repo_root <- resolve_repo_root()
pbmc_dir <- file.path(repo_root, "intermediate", "pbmc")
sim_dir <- file.path(pbmc_dir, "anndata_elements")
fig_dir <- file.path(repo_root, "output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

get_gene_names <- function(df) {
  rn <- rownames(df)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))) {
    return(as.character(rn))
  }
  
  candidate_cols <- c("gene", "genes", "gene_name", "gene_names", "symbol", "feature", "index", "X")
  hit <- intersect(candidate_cols, colnames(df))
  if (length(hit) > 0L) {
    return(as.character(df[[hit[1L]]]))
  }
  
  as.character(df[[1L]])
}

# singlet data
cts <- Matrix::readMM(file.path(pbmc_dir, "pbmc_myeloid_platelet_dbl_counts.mtx"))
var_df <- read.csv(file.path(pbmc_dir, "pbmc_myeloid_platelet_dbl_var.csv"), row.names = 1)
obs_df <- read.csv(file.path(pbmc_dir, "pbmc_myeloid_platelet_dbl_obs.csv"))

orig_gene_names <- get_gene_names(var_df)
if (length(orig_gene_names) != ncol(cts)) {
  stop("Original var_df gene names length does not match ncol(cts).")
}
colnames(cts) <- orig_gene_names
cts <- as(cts, "CsparseMatrix")

if (!("ann_types" %in% colnames(obs_df))) {
  stop("'ann_types' not found in obs_df.")
}

# sim data
sim_cts <- Matrix::readMM(file.path(sim_dir, "adata_pbmc_counts_mpa_sim.mtx"))
sim_var_df <- read.csv(file.path(sim_dir, "adata_pbmc_var_mpa_sim.csv"))
sim_obs_df <- read.csv(file.path(sim_dir, "adata_pbmc_obs_mpa_sim.csv"))

sim_gene_names <- get_gene_names(sim_var_df)
if (length(sim_gene_names) != ncol(sim_cts)) {
  stop("Simulated var_df gene names length does not match ncol(sim_cts).")
}
colnames(sim_cts) <- sim_gene_names
sim_cts <- as(sim_cts, "CsparseMatrix")

if (!("droplet_type" %in% colnames(sim_obs_df))) {
  stop("'droplet_type' not found in sim_obs_df.")
}

# keep only MPA_sim from sim data and rename simMPA
sim_keep <- !is.na(sim_obs_df$droplet_type) & sim_obs_df$droplet_type == "MPA_sim"
if (!any(sim_keep)) {
  stop("No rows found with droplet_type == 'MPA_sim'.")
}
sim_cts <- sim_cts[sim_keep, , drop = FALSE]
sim_obs_df <- sim_obs_df[sim_keep, , drop = FALSE]
sim_obs_df$ann_types <- "SimMPA"

# keep common genes; concatenate data
common_genes <- intersect(colnames(cts), colnames(sim_cts))
if (length(common_genes) < 100L) {
  stop("Too few intersecting genes: ", length(common_genes))
}

cts_common <- cts[, common_genes, drop = FALSE]
sim_cts_common <- sim_cts[, common_genes, drop = FALSE]
sim_cts_common <- sim_cts_common[, common_genes, drop = FALSE]

# ensure no cells in cts_common exist in sim_cts_common / tag line: duplicate, duplicated
stopifnot(sum(obs_df$barcode_2 %in% sim_obs_df$barcode_2)==0) # passes

cts_all <- rbind(cts_common, sim_cts_common)
cts_all <- as(cts_all, "CsparseMatrix")

# join meta; keep common columns
all_obs_cols <- union(colnames(obs_df), colnames(sim_obs_df))
for (nm in setdiff(all_obs_cols, colnames(obs_df))) obs_df[[nm]] <- NA
for (nm in setdiff(all_obs_cols, colnames(sim_obs_df))) sim_obs_df[[nm]] <- NA
obs_df <- obs_df[, all_obs_cols, drop = FALSE]
sim_obs_df <- sim_obs_df[, all_obs_cols, drop = FALSE]
obs_all <- rbind(obs_df, sim_obs_df)

# parameters
cluster_col <- "ann_types"
clusters_to_keep <- c("cDC1", "cDC2", "cMono", "Doublet", "MPA", "SimMPA", "nMono", "pDC", "Platelet")

# checks and subsetting
stopifnot(nrow(cts_all) == nrow(obs_all))
stopifnot(cluster_col %in% colnames(obs_all))

cluster_labels <- as.character(obs_all[[cluster_col]])
keep <- !is.na(cluster_labels) & cluster_labels %in% clusters_to_keep

cts_use <- cts_all[keep, , drop = FALSE]
obs_use <- obs_all[keep, , drop = FALSE]
cl_use <- droplevels(factor(cluster_labels[keep], levels = clusters_to_keep))

if (nlevels(cl_use) < 2L) stop("Need at least 2 clusters with cells.")

# aggregate
X <- Matrix::sparse.model.matrix(~ 0 + cl_use)   # cells x clusters
colnames(X) <- sub("^cl_use", "", colnames(X))
pb_counts <- Matrix::crossprod(cts_use, X)       # genes x clusters
print(Matrix::colSums(X))

# edgeR and voom
dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge, method = "TMM")

design0 <- matrix(1, ncol(pb_counts), 1)
# relax threshold for filtering
keep_genes <- filterByExpr(dge, design = design0, min.count = 1)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

design_voom <- matrix(1, ncol(dge), 1)
v <- voom(dge, design = design_voom, plot = FALSE)
E <- v$E  # genes x clusters

# whole transcriptome-informed correlation
cor_pearson <- cor(E, method = "pearson")

# dotplot; dendrogram calculated from cor_pearson
plot_cluster_dendrogram_fast <- function(
    cts_use,
    cell_meta,
    E,
    cor_pearson,
    cluster_col = "ann_types",
    gene_list = NULL,
    n_top_genes = 30,
    cluster_hclust_method = "complete",
    gene_hclust_method = "complete"
) {
  # dendrogram built on whole transcriptome; representative genes used for visualization
  if (!inherits(cts_use, "CsparseMatrix")) cts_use <- as(cts_use, "CsparseMatrix")
  if (!is.data.frame(cell_meta)) stop("cell_meta must be a data.frame.")
  if (!(cluster_col %in% colnames(cell_meta))) stop("cluster_col not found in cell_meta.")
  if (nrow(cts_use) != nrow(cell_meta)) stop("nrow(cts_use) must equal nrow(cell_meta).")
  if (is.null(colnames(cts_use))) stop("cts_use must have gene names in colnames.")
  if (is.null(rownames(E))) stop("E must have gene names in rownames.")
  
  # Whole-transcriptome dendrogram
  hc_clusters <- hclust(as.dist(1 - cor_pearson), method = cluster_hclust_method)
  cluster_order <- colnames(cor_pearson)[hc_clusters$order]
  cluster_order <- intersect(cluster_order, colnames(E))
  if (length(cluster_order) < 2L) stop("Need >=2 shared clusters between cor_pearson and E.")
  
  E_ordered <- E[, cluster_order, drop = FALSE]
  
  # Gene selection
  if (is.null(gene_list)) {
    gene_vars <- matrixStats::rowVars(E_ordered)
    ord <- order(gene_vars, decreasing = TRUE)
    n_take <- min(n_top_genes, length(ord))
    gene_selected <- rownames(E_ordered)[ord[seq_len(n_take)]]
    message("[plot_cluster_dendrogram_fast] Using top ", length(gene_selected), " variable genes.")
  } else {
    gene_selected <- intersect(as.character(gene_list), rownames(E_ordered))
    if (length(gene_selected) == 0L) stop("None of the provided genes were found in E rownames.")
    message("[plot_cluster_dendrogram_fast] Using ", length(gene_selected), " user-provided genes.")
  }
  
  # Optional gene ordering (within selected genes)
  E_sel <- E_ordered[gene_selected, , drop = FALSE]
  if (nrow(E_sel) >= 2L) {
    hc_genes <- hclust(
      as.dist(1 - cor(t(E_sel), use = "pairwise.complete.obs")),
      method = gene_hclust_method
    )
    gene_order <- rownames(E_sel)[hc_genes$order]
  } else {
    hc_genes <- NULL
    gene_order <- rownames(E_sel)
  }
  E_final <- E_sel[gene_order, , drop = FALSE]
  
  # Match genes to cells x genes matrix
  gene_idx <- match(gene_order, colnames(cts_use))
  ok <- !is.na(gene_idx)
  gene_idx <- gene_idx[ok]
  gene_order <- gene_order[ok]
  E_final <- E_final[gene_order, , drop = FALSE]
  if (length(gene_idx) == 0L) stop("Selected genes not found in cts_use colnames.")
  
  cts_sel <- cts_use[, gene_idx, drop = FALSE]
  
  # Cluster membership
  cl <- factor(as.character(cell_meta[[cluster_col]]), levels = cluster_order)
  keep_cells <- !is.na(cl)
  if (!all(keep_cells)) {
    cts_sel <- cts_sel[keep_cells, , drop = FALSE]
    cl <- droplevels(cl[keep_cells])
  }
  
  X <- Matrix::sparse.model.matrix(~ 0 + cl)
  colnames(X) <- sub("^cl", "", colnames(X))
  clusters_present <- colnames(X)
  cells_per_cluster <- Matrix::colSums(X)
  
  # Fraction expressing via sparse algebra
  nz <- cts_sel
  nz@x[] <- 1
  expr_counts <- Matrix::crossprod(nz, X)    # genes x clusters
  frac_expr <- sweep(expr_counts, 2, cells_per_cluster, "/", check.margin = FALSE)
  
  # Mean expression from voom matrix
  mean_log_cpm <- E_final[, clusters_present, drop = FALSE]
  
  # Long table
  df_frac <- as.data.frame(as.table(as.matrix(frac_expr)))
  colnames(df_frac) <- c("gene", "cluster", "frac_expressing")
  
  df_expr <- as.data.frame(as.table(as.matrix(mean_log_cpm)))
  colnames(df_expr) <- c("gene", "cluster", "mean_log_cpm")
  
  summary_df <- merge(df_frac, df_expr, by = c("gene", "cluster"), all = FALSE, sort = FALSE)
  
  summary_df$mean_expr_scaled <- ave(
    summary_df$mean_log_cpm,
    summary_df$gene,
    FUN = function(x) {
      rng <- range(x, na.rm = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] <= rng[1]) {
        rep(0, length(x))
      } else {
        (x - rng[1]) / (rng[2] - rng[1])
      }
    }
  )
  
  scanpy_reds_256 <- colorRampPalette(c(
    "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272",
    "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"
  ))(256)
  
  summary_df$gene <- factor(summary_df$gene, levels = gene_order, ordered = TRUE)
  summary_df$cluster <- factor(summary_df$cluster, levels = clusters_present, ordered = TRUE)
  summary_df$cluster_id <- as.numeric(summary_df$cluster)
  
  cluster_levels <- levels(summary_df$cluster)
  n_clust <- length(cluster_levels)
  
  # Pre-transform fraction with size_exponent = 1.5 (Scanpy default)
  summary_df$frac_size <- summary_df$frac_expressing ^ 1.5
  
  # Main dot plot
  p_main <- ggplot(
    summary_df,
    aes(x = gene, y = cluster_id, size = frac_size, fill = mean_expr_scaled)
  ) +
    geom_point(alpha = 0.85, pch = 21, color = "black", stroke = 0.2) + 
    scale_y_continuous(
      breaks = seq_len(n_clust),
      labels = cluster_levels,
      limits = c(0.5, n_clust + 0.5),
      expand = c(0, 0)
    ) + 
    scale_size_continuous(
      name   = "Fraction of cells\nin group (%)",
      range  = c(0, 9),              # 0 minimum = dots fade to nothing
      limits = c(0, 1),              # limits on the transformed scale
      breaks = c(0.2, 0.4, 0.6, 0.8, 1.0) ^ 1.5,   # transform legend breaks too
      labels = c(20, 40, 60, 80, 100)
    ) +
    scale_fill_gradientn(
      name   = "Mean expression\nin group",
      colors = scanpy_reds_256,
      limits = c(0, 1),
      guide  = guide_colorbar(
        frame.colour    = "black",
        frame.linewidth = 0.4,
        ticks.colour    = "black"
      )
    ) + 
    theme_bw(base_size = 20) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12), 
      legend.position = 'bottom'
    )
  
  # Right-side dendrogram (built from whole-transcriptome hc_clusters)
  dd <- ggdendro::dendro_data(as.dendrogram(hc_clusters), type = "rectangle")
  seg <- dd$segments
  
  # map dendrogram leaf positions (1..k in hc order) to plot y positions
  hc_labels_in_order <- hc_clusters$labels[hc_clusters$order]
  leaf_to_y <- setNames(match(hc_labels_in_order, cluster_levels), seq_along(hc_labels_in_order))
  
  map_leaf_pos <- function(v) {
    out <- v
    is_leaf <- abs(v - round(v)) < 1e-8
    out[is_leaf] <- leaf_to_y[as.character(round(v[is_leaf]))]
    out
  }
  
  seg$y_plot <- map_leaf_pos(seg$x)
  seg$yend_plot <- map_leaf_pos(seg$xend)
  
  p_dend <- ggplot(seg) +
    geom_segment(
      aes(x = y, y = y_plot, xend = yend, yend = yend_plot),
      linewidth = 0.6,
      lineend = "square"
    ) +
    scale_y_continuous(limits = c(0.5, n_clust + 0.5), expand = c(0, 0)) +
    scale_x_reverse(expand = c(0, 0)) +
    theme_void()
  
  p <- p_main + p_dend + patchwork::plot_layout(widths = c(12, 2))
  
  list(
    plot = p,
    summary_data = summary_df,
    cluster_order = clusters_present,
    gene_order = gene_order,
    cluster_dendrogram = hc_clusters,
    gene_dendrogram = hc_genes
  )
}

custom_genes <- c(
  "PPBP", "CAVIN2", "CLU", "SPARC", "GP9", "PRKAR2B", "HIST1H2AC",
  "ARHGAP18", "TRIM58", "S100A9", "S100A8", "LYZ", "CD14",
  "FCGR3A", "CLEC9A", "CLEC10A", "PLD4", "IRF4", "CD3E", "CD8A", "CD79A"
)

result_custom <- plot_cluster_dendrogram_fast(
  cts_use = cts_use,
  cell_meta = obs_use,
  E = E, 
  cor_pearson = cor_pearson,
  cluster_col = "ann_types",
  gene_list = custom_genes
)
# result_custom$plot
ggsave(filename = 'pbmc_myeloid_bubble.pdf', plot = result_custom$plot, device = cairo_pdf, 
  path = fig_dir, width = 9, height = 6)

