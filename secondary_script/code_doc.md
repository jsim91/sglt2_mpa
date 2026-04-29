# Code Documentation

This document describes the canonical execution sequence for scripts in `primary_script`, along with practical context for each step: role, key inputs, key outputs, dependencies, and handoff points.

## Start Here

- Workflow entrypoint: run `primary_script/1_score_solo_singlets.ipynb` first.
- That notebook begins the active pipeline by generating per-lane SOLO outputs consumed by later scripts.
- Typical first-pass run order is: `1` -> `2` -> `3` -> `4` -> `5` -> `6` -> `7` -> `8` -> `9_13` (step 9 pass) -> `10` -> `11` -> `12` -> `9_13` (step 13 pass) -> `14` -> `15`.
- Caveat: `9_13[FigA]_test_mast.R` is a combined script that supports two logical checkpoints in the run sequence (step 9 and step 13), so treat it as a re-entry script rather than a one-time step.

## Required H5 Dependency

- This workflow depends on external 10x Cell Ranger filtered matrix files (`filtered_feature_bc_matrix.h5`) for each lane.
- Expected location pattern in this repo: `primary_dependents/cellranger_h5/10872-MM-*/filtered_feature_bc_matrix.h5`.
- These large input files are not expected to be versioned in git and must be provisioned before running the entrypoint notebook.

## Data Availability Reference

- For data sharing policy and request language, see the Data Availability Statement in [README.md](../README.md#data-availability-statement).
- If raw or intermediate inputs are missing locally, use that statement as the source-of-truth for how to request access.

## Scope And Conventions

- This document covers only files in `primary_script`.
- Numeric prefixes indicate intended run order.
- `9_13[FigA]_test_mast.R` combines steps 9 and 13.
- Bracket tags (for example `[FigB]`) indicate manuscript figure linkage.
- Paths are relative to repository root (resolved by `REPO_ROOT` in notebooks and `here::here(...)` in R scripts).

## Pipeline Summary

1. Build per-lane singlet/doublet calls and merged PBMC objects (`1` to `3`).
2. Perform global PBMC integration and subset extraction (`4`).
3. Perform targeted myeloid integration and export downstream analysis tables (`5`).
4. Reintroduce and analyze doublet-inclusive populations (`6` and `7`).
5. Simulate MPA-like doublets and compare to observed MPAs (`8`).
6. Run downstream DGE/GSEA and figure-generation scripts (`9_13`, `10`, `11`, `12`, `14`).

## Script Sequence Details

### 1_score_solo_singlets.ipynb

- Role / Purpose:
	- Scores intra-sample doublets using SOLO (scvi-tools) on each lane after restricting to Souporcell singlets.
- Primary Inputs:
	- `primary_dependents/cellranger_h5/10872-MM-*/filtered_feature_bc_matrix.h5`
	- `primary_dependents/souporcell/10872-MM-*_soc/clusters.tsv`
- Primary Outputs:
	- `intermediate/solo_scores_MM1.csv` through `intermediate/solo_scores_MM8.csv`
- Other Relevant Context:
	- Loops through lanes 1 to 8.
	- Applies cell/gene filters (`min_genes=200`, `min_cells=10`) before model training.
	- These SOLO outputs are consumed by script 2 and script 6.

### 2_build_seurat.R

- Role / Purpose:
	- Builds lane-level Seurat objects, maps barcode metadata (Souporcell, SNP assignment, SOLO), applies QC, merges lanes, and exports matrix files.
- Primary Inputs:
	- Cell Ranger matrices in `primary_dependents/cellranger_h5`
	- Souporcell clusters in `primary_dependents/souporcell`
	- SOLO files from `intermediate/solo_scores_*`
	- SNP correlation files in `primary_dependents/snp_correlation`
	- Key file `primary_dependents/MM148_key.csv`
- Primary Outputs:
	- `intermediate/seurat/all_cells/mm_pbmc_lane_qc_metrics.rds`
	- `intermediate/seurat/all_cells/mm_pbmc_lane_list.rds`
	- `intermediate/seurat/all_cells/mm_pbmc_merged.rds`
	- `intermediate/seurat/RNA_assay_counts.mtx`
	- `intermediate/seurat/seurat_meta.csv`
	- `intermediate/seurat/gene_names.csv`
- Other Relevant Context:
	- Enforces Seurat and SeuratObject version checks.
	- Removes non-singlets and high-mito cells.
	- Exports the matrix/meta triplet expected by script 3.

### 3_convert_seurat_to_anndata.ipynb

- Role / Purpose:
	- Converts Seurat-exported matrix and metadata to an AnnData object.
- Primary Inputs:
	- `intermediate/seurat/RNA_assay_counts.mtx`
	- `intermediate/seurat/seurat_meta.csv`
	- `intermediate/seurat/gene_names.csv`
- Primary Outputs:
	- `intermediate/mm_pbmc.h5ad`
- Other Relevant Context:
	- This is the handoff into scanpy/scvi integration workflows.
	- Script 4 starts from this `h5ad` file.

### 4_integrate_global.ipynb

- Role / Purpose:
	- Runs global PBMC integration, clustering, broad cell-type annotation, and exports global PBMC artifacts.
- Primary Inputs:
	- `intermediate/mm_pbmc.h5ad`
	- `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
	- Cell metadata columns (lane, study_id, QC metrics)
- Primary Outputs:
	- `intermediate/pbmc/pbmc_int.h5ad`
	- `intermediate/pbmc/pbmc_final.h5ad`
	- `intermediate/pbmc/pbmc_subset_tnk.h5ad`
	- `intermediate/pbmc/pbmc_subset_b.h5ad`
	- `intermediate/pbmc/pbmc_subset_myeloid.h5ad`
	- `intermediate/pbmc/pbmc_subset_platelet.h5ad`
	- `intermediate/pbmc/pbmc_subset_hspc.h5ad`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_var.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_latent_coordinates.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_umap_coordinates.csv`
- Other Relevant Context:
	- Uses SCVI tuning/training, then neighbors/Leiden/UMAP.
	- Provides key inputs for scripts 5, 8, 9_13, and 12.

### 5_integrate_myeloid.ipynb

- Role / Purpose:
	- Re-integrates and refines the myeloid subset, then exports myeloid-specific annotations and matrices.
- Primary Inputs:
	- `intermediate/pbmc/pbmc_subset_myeloid.h5ad`
	- `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
- Primary Outputs:
	- `intermediate/pbmc/pbmc_myeloid_int.h5ad`
	- `intermediate/pbmc/pbmc_myeloid.h5ad`
	- `intermediate/pbmc/pbmc_myeloid_cell_map.csv`
	- `intermediate/pbmc/myeloid/object/pbmc_myeloid_small.h5ad`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_latent_coordinates.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv`
- Other Relevant Context:
	- Generates myeloid labels used downstream in scripts 7, 9_13, and bubble/figure workflows.

### 6_dbl_integrate_global.ipynb

- Role / Purpose:
	- Builds a doublet-inclusive global PBMC object (Souporcell and SOLO labels retained) and integrates/clusters it.
- Primary Inputs:
	- `primary_dependents/cellranger_h5/**`
	- `primary_dependents/souporcell/**`
	- `primary_dependents/snp_correlation/**`
	- `intermediate/solo_scores_MM*.csv`
	- `primary_dependents/MM148_key.csv`
- Primary Outputs:
	- `intermediate/pbmc/pbmc_int_dbl.h5ad`
	- `intermediate/pbmc/pbmc_myeloid_dbl.h5ad`
- Other Relevant Context:
	- Reconstructs lane-level metadata and doublet annotations directly from raw and intermediate sources.
	- Feeds script 7 for focused myeloid/platelet doublet analysis.

### 7_dbl_integrate_myeloid_platelet.ipynb

- Role / Purpose:
	- Focuses on myeloid + platelet populations from the doublet-inclusive object and computes marker contrasts.
- Primary Inputs:
	- `intermediate/pbmc/pbmc_int_dbl.h5ad`
	- `intermediate/pbmc/pbmc_final.h5ad` (for MPA singlet reference checks)
- Primary Outputs:
	- `intermediate/pbmc/pbmc_myeloid_platelet_int_dbl.h5ad`
	- `intermediate/pbmc/MPA_vs_Platelet_scanpy_markers.csv`
	- `intermediate/pbmc/MPA_vs_CD16Mono_scanpy_markers.csv`
	- `intermediate/pbmc/MPA_vs_CD14Mono_scanpy_markers.csv`
	- `intermediate/pbmc_myeloid_platelet_int_dbl_umap_coordinates_to_r.csv`
	- `intermediate/pbmc_myeloid_platelet_int_dbl_obs_to_r.csv`
- Other Relevant Context:
	- Exports R-friendly files later consumed by script 9_13 for figure panels.
	- Establishes cluster-level doublet frequency metrics for comparison to MPA cluster.

### 8_dbl_simulate_myeloid_platelet.ipynb

- Role / Purpose:
	- Simulates cMono+platelet doublets (SimMPA), integrates with observed cells, and exports simulation artifacts.
- Primary Inputs:
	- `intermediate/pbmc/pbmc_myeloid.h5ad`
	- `intermediate/pbmc/pbmc_final.h5ad`
	- `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
- Primary Outputs:
	- `intermediate/pbmc/pbmc_simulation_dataset.h5ad`
	- `intermediate/pbmc/pbmc_myeloid_platelet_int_sim.h5ad`
	- `intermediate/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_counts_mpa_sim.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_var_mpa_sim.csv`
- Other Relevant Context:
	- SimMPA cells are tagged via barcode suffix logic.
	- Outputs are consumed by scripts 9_13 and 11.

### 9_13[FigA]_test_mast.R (Steps 9 and 13; Figure A)

- Role / Purpose:
	- Runs MAST DGE workflows for key PBMC/myeloid contrasts, creates panel-relevant visualizations, and writes intermediate MAST products.
- Primary Inputs:
	- `intermediate/pbmc/anndata_elements/adata_pbmc_*.{mtx,csv}`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_*.csv`
	- `intermediate/pbmc_myeloid_platelet_int_dbl_*_to_r.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv`
	- `intermediate/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates.csv`
	- `primary_dependents/seurat_dge_local.R`
- Primary Outputs:
	- `intermediate/seurat/MM148_pbmc_seurat.rds`
	- `intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds`
	- `intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_gsea_sd1sd3_15JAN2025.rds`
	- `intermediate/mast/platelet_genes/mast_platelet_gene_result.csv`
	- `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.csv`
	- `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.csv`
	- `intermediate/mast/mpa_sd_compare/mast_gsea_mpa_vs_cmono_sd1_29Jan2024.rds`
	- `intermediate/mast/mpa_sd_compare/mast_gsea_mpa_vs_cmono_sd3_29Jan2024.rds`
	- figure exports in `output/figures`
- Other Relevant Context:
	- This script provides critical inputs used directly by scripts 11 and 14.
	- Includes platelet-gene blacklist derivation for cleaner MPA versus cMono contrasts.

### 10[FigB]_plot_bubble.ipynb (Figure B)

- Role / Purpose:
	- Builds the marker bubble visualization panel across monocyte, platelet, MPA, simulated, and doublet populations.
- Primary Inputs:
	- `intermediate/pbmc/pbmc_final.h5ad`
	- `intermediate/pbmc/pbmc_myeloid.h5ad`
	- `intermediate/pbmc/pbmc_myeloid_platelet_int_dbl.h5ad`
	- simulation object(s) from script 8
- Primary Outputs:
	- Bubble plot figure exports (for manuscript panel B)
- Other Relevant Context:
	- This script is figure-focused and typically run after core object generation is complete.

### 11_test_mast_real_vs_sim.R (Upstream Input For Script 14)

- Role / Purpose:
	- Performs MAST DGE for real MPA versus simulated MPA and applies platelet-gene blacklist filtering.
- Primary Inputs:
	- `intermediate/pbmc/anndata_elements/adata_pbmc_counts_mpa_sim.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_var_mpa_sim.csv`
	- `intermediate/mast/platelet_genes/mast_platelet_gene_result.csv`
	- `primary_dependents/seurat_dge_local.R`
- Primary Outputs:
	- `intermediate/mast/mpa_real_sim/pbmc_mpa_sim_real_no_platelet_genes_deg_6FEB2025.rds`
- Other Relevant Context:
	- This output is a direct input to script 14.

### 12_summarize_anndata_to_seurat.R

- Role / Purpose:
	- Reconstructs Seurat object from exported AnnData components and creates sample-by-cluster count summaries for downstream visualization scripts.
- Primary Inputs:
	- `intermediate/pbmc/anndata_elements/adata_pbmc_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_var.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_umap_coordinates.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_latent_coordinates.csv`
- Primary Outputs:
	- `intermediate/pbmc/anndata_elements/pbmc_seurat.rds`
	- `intermediate/pbmc/anndata_elements/mm_pbmc_subset_cluster_cell_counts.csv`
	- `intermediate/pbmc/anndata_elements/mm_pbmc_named_cluster_cell_counts.csv`
- Other Relevant Context:
	- These count tables are used for downstream composition and correlation analyses.

### 14[FigC_H_I]_run_fgsea.R (Figure C, H, I)

- Role / Purpose:
	- Runs fgsea on ranked MAST contrasts and writes figure-ready enrichment plots.
- Primary Inputs:
	- `intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds`
	- `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.csv`
	- `intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.csv`
	- `intermediate/mast/mpa_real_sim/pbmc_mpa_sim_real_no_platelet_genes_deg_6FEB2025.rds`
	- `primary_dependents/gene_set/H_Hallmark_modules_gsea.gmt`
- Primary Outputs:
	- `output/figures/sim_dblt_wo_plt_gsea.pdf`
	- `output/figures/combined_sd1_sd3_gsea.pdf`
	- `output/figures/mpa_cmono_gsea.pdf`
- Other Relevant Context:
	- Depends on outputs from scripts 9_13 and 11.
	- Uses Hallmark pathways and NES-based visualization.

### 15[FigE]_circ_mm148.R (Figure E)

- Role / Purpose:
	- Analyzes and visualizes cluster frequencies (MPA and cMono) across study timepoints using violin plots and publication-ready graphics.
- Primary Inputs:
	- `intermediate/pbmc/anndata_elements/mm_pbmc_named_cluster_cell_counts.csv` (from script 12)
- Primary Outputs:
	- Violin plot figures for MPA and cMono frequency visualization (displayed/saved during script execution)
- Other Relevant Context:
	- Consumes the count summary table generated by script 12.
	- Depends on script 12 completion before running.
	- Figure-generation script typically run after upstream data processing finalized.

## Dependency Handoff Map

- `1` produces SOLO scores used by `2` and `6`.
- `2` produces matrix/meta exports used by `3`.
- `3` produces global AnnData used by `4`.
- `4` produces PBMC subsets and anndata elements used by `5`, `8`, `9_13`, and `12`.
- `5` produces myeloid exports used by `7` and `9_13`.
- `6` produces doublet-inclusive PBMC object used by `7`.
- `7` produces myeloid/platelet doublet exports used by `9_13`.
- `8` produces simulation outputs used by `9_13` and `11`.
- `9_13` produces MAST outputs used by `11` and `14`.
- `11` produces real-versus-sim MAST output used by `14`.
- `12` produces count summaries used by `15`.

## Execution Notes

- Notebook scripts rely on a `REPO_ROOT` resolver and generally run from the repository tree.
- R scripts rely on `here::here(...)` and should be run with the repo root as project root.
- Some scripts include optional or cached read/write blocks. Verify whether to regenerate or reuse existing artifacts before reruns.
- Figure-tagged scripts (`10`, `11`, `14`, `15`, and portions of `9_13`) are usually run after upstream objects are finalized.
