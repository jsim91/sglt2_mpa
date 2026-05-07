# Code Documentation

This document describes the canonical execution sequence for scripts in `primary_script`, with emphasis on what is needed to reproduce the primary manuscript outputs.

## Start Here

- Workflow entrypoint: run `primary_script/1_pre_detect_doublets.ipynb` first.
- Primary sequence (single pass): `1` -> `2` -> `3` -> `4` -> `5` -> `6` -> `7` -> `8` -> `9` -> `10` -> `11` -> `12` -> `13` -> `14`.

## Before Running

- External 10x Cell Ranger H5 matrices are required and are not packaged in this repo.
- Placeholder directory structure is included, but files must be requested/provisioned before execution.
- Expected H5 pattern: `primary_dependents/cellranger_h5/10872-MM-*/filtered_feature_bc_matrix.h5`.
- Data sharing/request policy is in [README.md](../README.md#data-availability-statement).

## How To Run

- Notebooks: run from repository root so `REPO_ROOT` resolvers and relative paths work.
- R scripts: run from repository root (or an R project rooted at repo root) so `here::here(...)` resolves correctly.
- Suggested invocation style for R scripts:
  - `Rscript primary_script/2_pre_build_seurat.R`
  - `Rscript primary_script/9[FigA]_umaps_mast.R`

## Scope And Conventions

- Prefixes encode branch intent:
  - `pre` = shared preparation
  - `sgl` = singlet branch
  - `dbl` = doublet-inclusive branch
  - `sim` = simulation branch
- Scripts `9+` are mostly downstream/branch-agnostic analysis and plotting.
- Figure tags in filenames (`[FigB]`, `[FigC_H_I]`, etc.) indicate manuscript panel linkage.

## Script Sequence Details

### 1_pre_detect_doublets.ipynb

- Purpose:
  - Run SOLO intra-sample doublet scoring on Souporcell-singlet-filtered lanes.
- Primary inputs:
  - `primary_dependents/cellranger_h5/10872-MM-*/filtered_feature_bc_matrix.h5`
  - `primary_dependents/souporcell/10872-MM-*_soc/clusters.tsv`
- Primary outputs:
  - `intermediate/solo_scores_MM1.csv` ... `intermediate/solo_scores_MM8.csv`
- Context:
  - Starts the end-to-end pipeline by generating barcode-level doublet scores used in prep and doublet-inclusive branches.

### 2_pre_build_seurat.R

- Purpose:
  - Build lane-level and merged Seurat objects, attach genotype/doublet metadata, apply QC, and export matrix/meta triplet.
- Primary inputs:
  - Cell Ranger matrices (`primary_dependents/cellranger_h5`)
  - Souporcell clusters (`primary_dependents/souporcell`)
  - SOLO scores (`intermediate/solo_scores_*`)
  - SNP correlation files (`primary_dependents/snp_correlation`)
  - `primary_dependents/MM148_key.csv`
- Primary outputs:
  - `intermediate/seurat/RNA_assay_counts.mtx`
  - `intermediate/seurat/seurat_meta.csv`
  - `intermediate/seurat/gene_names.csv`
- Context:
  - Standardized R-side object construction and QC handoff into AnnData/scverse workflows.

### 3_pre_seurat_to_anndata.ipynb

- Purpose:
  - Convert Seurat exports to canonical PBMC AnnData.
- Primary inputs:
  - `intermediate/seurat/RNA_assay_counts.mtx`
  - `intermediate/seurat/seurat_meta.csv`
  - `intermediate/seurat/gene_names.csv`
- Primary output:
  - `intermediate/mm_pbmc.h5ad`
- Context:
  - Bridge from Seurat preprocessing into integration/clustering notebooks.

### 4_sgl_integrate_global.ipynb

- Purpose:
  - Global singlet-branch PBMC integration/clustering and broad lineage partitioning.
- Primary inputs:
  - `intermediate/mm_pbmc.h5ad`
  - `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
- Primary outputs (used downstream):
  - `intermediate/pbmc/pbmc_int.h5ad`
  - `intermediate/pbmc/pbmc_subset_myeloid.h5ad`
  - `intermediate/pbmc/pbmc_subset_platelet.h5ad`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_counts.mtx`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_obs.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_var.csv`
- Context:
  - Establishes the main singlet integrated reference and subset objects used by both doublet and simulation branches.

### 5_sgl_integrate_myeloid.ipynb

- Purpose:
  - Refine/integrate the singlet myeloid subset and export myeloid elements for R analyses.
- Primary inputs:
  - `intermediate/pbmc/pbmc_subset_myeloid.h5ad`
  - `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
- Primary outputs (used downstream):
  - `intermediate/pbmc/pbmc_myeloid.h5ad`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_latent_coordinates.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv`
- Context:
  - Defines canonical myeloid labels including MPA in the singlet branch.

### 6_dbl_integrate_global.ipynb

- Purpose:
  - Rebuild global PBMC with Souporcell doublets added back and re-integrate in the doublet-inclusive branch.
- Primary inputs:
  - Raw/metadata dependencies from `primary_dependents` (Cell Ranger, Souporcell, SNP, key)
  - `intermediate/pbmc/pbmc_int.h5ad`
  - `intermediate/pbmc/pbmc_myeloid.h5ad`
  - `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
- Primary outputs:
  - `intermediate/pbmc/pbmc_int_dbl.h5ad`
  - `intermediate/pbmc/pbmc_myeloid_dbl.h5ad`
- Context:
  - Creates the technical-doublet-aware branch used to test whether MPA behaves like artifact.

### 7_dbl_integrate_myeloid_platelet.ipynb

- Purpose:
  - Focused myeloid+platelet integration in the doublet-inclusive branch; compute and export doublet metrics.
- Primary inputs:
  - `intermediate/pbmc/pbmc_int_dbl.h5ad`
  - `intermediate/pbmc/pbmc_myeloid_dbl.h5ad`
  - `intermediate/pbmc/pbmc_myeloid.h5ad`
  - `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
- Primary outputs:
  - `intermediate/pbmc/pbmc_myeloid_platelet_int_dbl.h5ad`
  - `intermediate/pbmc/pbmc_myeloid_platelet_dbl_counts.mtx`
  - `intermediate/pbmc/pbmc_myeloid_platelet_dbl_obs.csv`
  - `intermediate/pbmc/pbmc_myeloid_platelet_dbl_var.csv`
  - `intermediate/pbmc_myeloid_platelet_int_dbl_umap_coordinates_to_r.csv`
  - `intermediate/pbmc_myeloid_platelet_int_dbl_obs_to_r.csv`
- Context:
  - Supplies doublet-branch matrices/metadata to Figure A and Figure B analysis scripts.

### 8_sim_simulate_MPA.ipynb

- Purpose:
  - Build simulated MPA proxy cells (SimMPA), integrate with observed cells, and export simulation artifacts.
- Primary inputs:
  - `intermediate/pbmc/pbmc_myeloid.h5ad`
  - `intermediate/pbmc/pbmc_subset_platelet.h5ad`
  - `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
- Primary outputs:
  - `intermediate/pbmc/pbmc_myeloid_platelet_int_sim.h5ad`
  - `intermediate/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_counts_mpa_sim.mtx`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_var_mpa_sim.csv`
- Context:
  - Simulation branch to benchmark observed MPAs against controlled SimMPA profiles.

### 9[FigA]_umaps_mast.R (Figure A + MAST Hand-offs)

- Purpose:
  - Generate three Figure A UMAP panels (singlet, doublet-inclusive, simulation) and run key MAST contrasts.
- Primary inputs:
  - Myeloid/platelet anndata elements from scripts 4 and 5
  - Doublet-branch exports from script 7
  - Simulation exports from script 8
  - `primary_dependents/seurat_dge_local.R`
- Primary outputs (used downstream):
  - `intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds`
  - `intermediate/mast/platelet_genes/mast_platelet_gene_result.csv`
  - `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.csv`
  - `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.csv`
- Additional outputs (figure artifacts):
  - `output/figures/pbmc_myeloid_singlet_umap_16APR2025.pdf`
  - `output/figures/pbmc_myeloid_doublet_umap_16APR2025.pdf`
  - `output/figures/pbmc_myeloid_sim_umap_16APR2025.pdf`
- Context:
  - Central handoff script for MAST inputs consumed by scripts 11 and 13.

### 10[FigB]_bulk_gene_dotplot.R (Figure B)

- Purpose:
  - Build the Figure B bulk-gene dotplot comparing observed myeloid/platelet states and SimMPA context.
- Primary inputs:
  - `intermediate/pbmc/pbmc_myeloid_platelet_dbl_counts.mtx`
  - `intermediate/pbmc/pbmc_myeloid_platelet_dbl_obs.csv`
  - `intermediate/pbmc/pbmc_myeloid_platelet_dbl_var.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_counts_mpa_sim.mtx`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_var_mpa_sim.csv`
- Primary output:
  - `output/figures/pbmc_myeloid_bubble.pdf`
- Context:
  - Figure B visualization script.

### 11_test_mast_real_vs_sim.R

- Purpose:
  - Run MAST DGE contrast of `MPA_real` vs `MPA_sim` with platelet-gene blacklist filtering.
- Primary inputs:
  - Simulation anndata elements from script 8
  - `intermediate/mast/platelet_genes/mast_platelet_gene_result.csv` from script 9
  - `primary_dependents/seurat_dge_local.R`
- Primary output:
  - `intermediate/mast/mpa_real_sim/pbmc_mpa_sim_real_no_platelet_genes_deg_6FEB2025.rds`
- Context:
  - Produces the real-vs-sim ranked contrast used by script 13.

### 12_write_celltype_counts.R

- Purpose:
  - Write myeloid cell-type count table by participant/timepoint for frequency plotting.
- Primary inputs:
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv`
  - `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_latent_coordinates.csv`
- Primary output (downstream-relevant):
  - `intermediate/pbmc/anndata_elements/mm_myeloid_named_cluster_cell_counts.csv`
- Context:
  - Supports Figure E frequency plotting in script 14.

### 13[FigC_H_I]_plot_fgsea.R (Figure C, H, I)

- Purpose:
  - Plot fgsea outputs for SD comparisons and real-vs-sim contrasts shown in panels C, H, and I.
- Primary inputs:
  - `intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds`
  - `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.csv`
  - `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.csv`
  - `intermediate/mast/mpa_real_sim/pbmc_mpa_sim_real_no_platelet_genes_deg_6FEB2025.rds`
  - `primary_dependents/gene_set/H_Hallmark_modules_gsea.gmt`
- Primary outputs:
  - `output/figures/combined_sd1_sd3_gsea.pdf`
  - `output/figures/mpa_cmono_gsea.pdf`
  - `output/figures/sim_dblt_wo_plt_gsea.pdf`
- Context:
  - Converts ranked DGE contrasts into figure-ready pathway enrichment panels.

### 14[FigE]_plot_myeloid_frequency.R (Figure E)

- Purpose:
  - Plot myeloid frequency trajectories (notably MPA and cMono) across study days.
- Primary input:
  - `intermediate/pbmc/anndata_elements/mm_myeloid_named_cluster_cell_counts.csv`
- Primary output:
  - `output/figures/mpa_cmono_frequency.pdf`
- Context:
  - Final composition-style visualization panel linked to manuscript Figure E.

## Dependency Handoff (Primary Sequence)

- `1` -> `2`: SOLO scores.
- `2` -> `3`: matrix/meta/gene exports.
- `3` -> `4`: global PBMC AnnData.
- `4` -> `5`: myeloid subset; `4` -> `8` platelet subset; `4` -> `9` platelet anndata elements.
- `5` -> `6,7,8,9,12`: myeloid object and anndata elements.
- `6` -> `7`: doublet-inclusive PBMC/myeloid objects.
- `7` -> `9,10`: doublet-branch matrix/metadata exports.
- `8` -> `9,10,11`: simulation matrix/metadata exports.
- `9` -> `11,13`: MAST outputs and platelet blacklist.
- `11` -> `13`: real-vs-sim MAST output.
- `12` -> `14`: cell-type count table.

## Practical Notes

- Run notebooks and R scripts from repo root for stable path resolution.
- `9[FigA]_umaps_mast.R` uses an automatic cache for `intermediate/seurat/MM148_pbmc_seurat.rds` (build if missing, read if present).
