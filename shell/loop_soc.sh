#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BARCODE_SCRIPT="${REPO_ROOT}/secondary_script/retrieve_barcodes.R"
CELLRANGER_ROOT="${REPO_ROOT}/primary_dependents/cellranger_h5"
SOUPORCELL_OUT_ROOT="${REPO_ROOT}/primary_dependents/souporcell"

# Override with env vars if your .sif and reference fasta are elsewhere.
SINGULARITY_BIN="${SINGULARITY_BIN:-singularity}"
SOUPORCELL_SIF="${SOUPORCELL_SIF:-${REPO_ROOT}/primary_dependents/souporcell_latest.sif}"
REFERENCE_FASTA="${REFERENCE_FASTA:-${REPO_ROOT}/primary_dependents/hg38.fa}"

for LANE in {1..8}
do
    echo "Processing lane ${LANE}"

    LANE_DIR="${CELLRANGER_ROOT}/10872-MM-${LANE}"
    OUT_DIR="${SOUPORCELL_OUT_ROOT}/10872-MM-${LANE}_soc"

    Rscript "${BARCODE_SCRIPT}" "${LANE_DIR}"

    nohup "${SINGULARITY_BIN}" exec -B "${LANE_DIR}" "${SOUPORCELL_SIF}" \
        souporcell_pipeline.py \
        -i "${LANE_DIR}/possorted_genome_bam.bam" \
        -b "${LANE_DIR}/barcodes_R.tsv" \
        -f "${REFERENCE_FASTA}" \
        -t 40 \
        -o "${OUT_DIR}" \
        -k 3 &
done

exit 0
