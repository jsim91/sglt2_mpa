Per-lane count matrices extracted from 10x Cell Ranger HDF5 files.

Expected files (one set per lane, lanes MM-1 through MM-8):
  10872-MM-{1..8}_counts.mtx     -- sparse count matrix (genes x barcodes)
  10872-MM-{1..8}_obs_names.csv  -- barcode (cell) names
  10872-MM-{1..8}_var_names.csv  -- feature (gene) names

These files are consumed by primary_script/7_score_doublets_scds.R.

These files are not versioned in git and must be generated before running
script 7.  They can be produced from the Cell Ranger HDF5 files located in
primary_dependents/cellranger_h5/ using the provided secondary_script utility,
or requested from the study group (see README.md Data Availability Statement).
