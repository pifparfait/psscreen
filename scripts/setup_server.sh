#!/usr/bin/env bash
# One-shot setup for a Linux server (no root needed): installs Miniforge if conda is missing,
# creates the `psscreen` environment, and runs the quick tests.
# Usage: bash scripts/setup_server.sh
set -euo pipefail

if ! command -v conda >/dev/null 2>&1; then
  echo ">> conda not found: installing Miniforge into \$HOME/miniforge3"
  url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  curl -L -o /tmp/miniforge.sh "$url"
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
  # shellcheck disable=SC1091
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
else
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi

cd "$(dirname "$0")/.."
conda env create -f environment.yml -n psscreen || conda env update -f environment.yml -n psscreen
conda activate psscreen
echo ">> machine: $(nproc) cores, $(free -g | awk '/Mem:/{print $2}') GB RAM"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
python tests/test_soc.py && python tests/test_tdwin.py
echo ">> OK. Next: bash scripts/run_validation_server.sh"
