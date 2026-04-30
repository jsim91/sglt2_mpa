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
	- Runs global PBMC integration, clustering, broad cell-type annotation, exports cell-type subsets, and writes platelet matrix elements for downstream R analysis.
- Primary Inputs:
	- `intermediate/mm_pbmc.h5ad`
	- `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
	- Cell metadata columns (lane, study_id, QC metrics)
- Primary Outputs:
	- `intermediate/pbmc/pbmc_subset_tnk.h5ad`
	- `intermediate/pbmc/pbmc_subset_b.h5ad`
	- `intermediate/pbmc/pbmc_subset_myeloid.h5ad`
	- `intermediate/pbmc/pbmc_subset_platelet.h5ad`
	- `intermediate/pbmc/pbmc_subset_hspc.h5ad`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_var.csv`
- Other Relevant Context:
	- Uses SCVI tuning/training, then neighbors/Leiden/UMAP.
	- Provides key inputs for scripts 5, 8, 9_13, and 10.
- Cell Annotation:
	- Leiden clustering at resolution=3 produces integer overclusters.
	- Overclusters are manually mapped to broad lineage labels via `merge_type_ref` dictionary (inspected via marker dotplot): TNK (clusters 0,2–4,7,9–12,14–18,22,25,27,29,35–37), B (1,13,20,30,31,33), Myeloid (5,6,8,19,21,23,26,28,34), Platelet (24), HSPC (32).
	- Label stored in `merged_type` column; used only to partition subsets for downstream scripts.

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
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_latent_coordinates.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv`
- Other Relevant Context:
	- Generates myeloid labels used downstream in scripts 6, 9_13, and bubble/figure workflows.
- Cell Annotation:
	- Leiden clustering at resolution=1 produces integer overclusters on the myeloid subset only.
	- Overclusters are manually mapped to refined myeloid labels via `myeloid_types` dictionary (inspected via marker dotplot): CD14 Mono (1–6,8,10–13), CD16 Mono (0), M-platelet (7), cDC2 (9), pDC (14), cDC1 (15).
	- Label stored in `merged_type` column; this is the canonical myeloid cell-type identity used by all downstream scripts.

### 6_dbl_integrate_global.ipynb

- Role / Purpose:
	- Builds a doublet-inclusive global PBMC object (Souporcell and SOLO labels retained) and integrates/clusters it.
- Primary Inputs:
	- `primary_dependents/cellranger_h5/**`
	- `primary_dependents/souporcell/**`
	- `primary_dependents/snp_correlation/**`
	- `intermediate/solo_scores_MM*.csv`
	- `primary_dependents/MM148_key.csv`
	- `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
	- `intermediate/pbmc/pbmc_myeloid.h5ad` (for M-platelet barcode identification)
- Primary Outputs:
	- `intermediate/pbmc/pbmc_int_dbl.h5ad`
	- `intermediate/pbmc/pbmc_myeloid_dbl.h5ad`
- Other Relevant Context:
	- Reconstructs lane-level metadata and doublet annotations directly from raw and intermediate sources.
	- Uses the singlet myeloid object (Script 5) to identify M-platelet barcodes.
	- Feeds script 7 for focused myeloid/platelet doublet analysis.
- Cell Annotation:
	- Leiden clustering at resolution=3 is run on all cells including Souporcell and SOLO doublets, producing integer `overcluster` labels.
	- Overclusters are mapped to broad lineage groups via a `type_ref` dictionary; Myeloid clusters (6, 9, 12, 16, 17, 21, 23, 24, 26, 27, 28, 33, 34, 40) and Platelet cluster (31) are identified by marker inspection. Label is stored in `cell_group`.
	- M-platelet barcodes from Script 5 (`merged_type == 'M-platelet'`) are used to verify which overclusters in the doublet-inclusive space contain the singlet MPA population.
	- The exported `pbmc_myeloid_dbl.h5ad` contains only cells where `cell_group == 'Myeloid'`, capturing both true myeloid singlets and doublets that co-cluster with myeloid populations.

### 7_dbl_integrate_myeloid_platelet.ipynb

- Role / Purpose:
	- Focuses on myeloid + platelet populations from the doublet-inclusive object and computes marker contrasts.
- Primary Inputs:
	- `intermediate/pbmc/pbmc_int_dbl.h5ad`
	- `intermediate/pbmc/pbmc_myeloid.h5ad` (singlet myeloid reference)
	- `primary_dependents/EXCLUDE_XY_TCR_IG.csv`
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
- Cell Annotation:
	- Leiden overclustering is run independently on the doublet-inclusive myeloid+platelet object (no labels carried from Script 5).
	- Overclusters are manually assigned to `ann_types` by inspecting marker expression: cMo (clusters 9,2,4,12,8,15,17,20,3,6), MPA (11), nMo (0), Platelet (16), pDC (18), cDC2 (13), cDC1 (21), doublet (19,1,22,5,14,7,10).
	- Per-cluster Souporcell and SOLO doublet rates are computed and stored as per-cell frequency metadata columns (`doublet_frequency`, `souporcell_doublet_frequency`, `doublet_sum_frequency`).
	- Script 5 `merged_type` is used only to look up which overclusters the known M-platelet barcodes fall in; it does not propagate as the working cell identity label.

### 8_dbl_simulate_myeloid_platelet.ipynb

- Role / Purpose:
	- Simulates cMono+platelet doublets (SimMPA), integrates with observed cells, and exports simulation artifacts.
- Primary Inputs:
	- `intermediate/pbmc/pbmc_myeloid.h5ad`
	- `intermediate/pbmc/pbmc_subset_platelet.h5ad`
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
- Cell Annotation:
	- Synthetic doublet cells are labelled `solo_prediction = 'cMo_Platelet_doublet'` and `souporcell_status = 'sim_doublet'` at construction time.
	- After scVI integration, Leiden clustering at resolution=2 produces new integer overclusters on the combined real+simulated object (independent of Scripts 5 and 7).
	- Overclusters are manually mapped to `cell_type` via `cell_type_dict`: CD14 Mono (2–12,15–19,21,22), CD16 Mono (0,11), MPA (1,26), Int Mono (3), Platelet (13), cDC2 (14), pDC (20), cDC1 (25), migratory myeloid (23), other myeloid (24).
	- A final `droplet_type` label is derived: cells with `cell_type == 'MPA'` → `MPA_real`; cells with `solo_prediction == 'cMo_Platelet_doublet'` → `MPA_sim`; all others inherit their `solo_prediction` value. This is the label used by scripts 9_13 and 11 for DGE.

### 9_13[FigA]_test_mast.R (Steps 9 and 13; Figure A)

- Role / Purpose:
	- Runs MAST DGE workflows for key PBMC/myeloid contrasts, creates panel-relevant visualizations, and writes intermediate MAST products.
- Primary Inputs:
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_platelet_var.csv`
	- `intermediate/pbmc_myeloid_platelet_int_dbl_obs_to_r.csv`
	- `intermediate/pbmc_myeloid_platelet_int_dbl_umap_coordinates_to_r.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv`
	- `intermediate/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates.csv`
	- `intermediate/seurat/MM148_pbmc_seurat.rds` (cached; rebuilt from elements when `if(F)` block is active)
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
	- `intermediate/pbmc/pbmc_subset_platelet.h5ad`
	- `intermediate/pbmc/pbmc_myeloid.h5ad`
	- `intermediate/pbmc/pbmc_myeloid_platelet_int_dbl.h5ad`
	- `intermediate/pbmc/pbmc_myeloid_platelet_int_sim.h5ad`
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
	- Reconstructs Seurat object from myeloid-refined AnnData components (Script 5 outputs) and creates myeloid cell-type-by-sample count summaries for downstream visualization.
- Primary Inputs:
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv`
	- `intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_latent_coordinates.csv`
- Primary Outputs:
	- `intermediate/pbmc/anndata_elements/pbmc_seurat.rds`
	- `intermediate/pbmc/anndata_elements/mm_myeloid_named_cluster_cell_counts.csv`
- Other Relevant Context:
	- Uses myeloid-specific refined annotations from Script 5 (merged_type: M-platelet, CD14 Mono, CD16 Mono, cDC1, cDC2, pDC).
	- Myeloid count table consumed by Script 15 for frequency visualization.

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

### 15[FigE]_plot_myeloid_frequency.R (Figure E)

- Role / Purpose:
	- Visualizes myeloid cell-type frequencies (MPA and cMono) across study timepoints using refined myeloid clustering from Script 5.
- Primary Inputs:
	- `intermediate/pbmc/anndata_elements/mm_myeloid_named_cluster_cell_counts.csv` (from Script 12)
- Primary Outputs:
	- `output/figures/mpa_cmono_frequency.pdf`
- Other Relevant Context:
	- Uses myeloid-normalized frequencies (% of total myeloid cells per sample).
	- Focuses on key comparisons: SD1→SD2 and SD1→SD3.
	- Figure-generation script typically run after upstream data processing finalized.

## Dependency Handoff Map

- `1` produces SOLO scores used by `2` and `6`.
- `2` produces matrix/meta exports used by `3`.
- `3` produces global AnnData used by `4`.
- `4` produces PBMC subsets and platelet anndata elements used by `5`, `8`, `9_13`, and `10`.
- `5` produces myeloid exports used by `6`, `7`, `9_13`, `12`, and `15`.
- `6` produces doublet-inclusive PBMC object used by `7`.
- `7` produces myeloid/platelet doublet exports used by `9_13` and `10`.
- `8` produces simulation outputs used by `9_13` and `11`.
- `9_13` produces MAST outputs used by `11` and `14`.
- `11` produces real-versus-sim MAST output used by `14`.
- `12` produces myeloid count summaries used by `15`.

## Execution Notes

- Notebook scripts rely on a `REPO_ROOT` resolver and generally run from the repository tree.
- R scripts rely on `here::here(...)` and should be run with the repo root as project root.
- Some scripts include optional or cached read/write blocks. Verify whether to regenerate or reuse existing artifacts before reruns.
- Figure-tagged scripts (`10`, `11`, `14`, `15`, and portions of `9_13`) are usually run after upstream objects are finalized.
