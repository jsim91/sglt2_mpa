#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_FILE="${REPO_ROOT}/environment.yml"
POST_INSTALL_R="${SCRIPT_DIR}/post_install_r_packages.R"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda is required but was not found on PATH" >&2
  exit 1
fi

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Could not find ${ENV_FILE}" >&2
  exit 1
fi

ENV_NAME="${1:-$(awk '/^name:/ {print $2; exit}' "${ENV_FILE}")}"
if [[ -z "${ENV_NAME}" ]]; then
  echo "Could not determine conda environment name" >&2
  exit 1
fi

echo "Repository root: ${REPO_ROOT}"
echo "Environment file: ${ENV_FILE}"
echo "Environment name: ${ENV_NAME}"

if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
  echo "Updating existing conda environment"
  conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
  echo "Creating conda environment"
  conda env create -f "${ENV_FILE}"
fi

echo "Installing SeuratWrappers from bu_cnio"
conda install -n "${ENV_NAME}" -y -c bu_cnio r-seuratwrappers

echo "Installing supplemental R packages"
conda run -n "${ENV_NAME}" Rscript "${POST_INSTALL_R}" "${ENV_NAME}"

echo "Environment setup complete"