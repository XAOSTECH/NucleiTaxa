#!/bin/bash
set -e

# post-create.sh
# Optional runtime setup: installs R packages and (optionally) Python packages.
# To skip Python installs during development/rebuild, set SKIP_PYTHON=1 in your environment.

if [[ "${SKIP_PYTHON:-0}" == "1" ]]; then
  echo "SKIP_PYTHON=1 - skipping Python package installs"
else
  if command -v pip3 >/dev/null 2>&1 && [[ -f "/workspaces/GAPMananas/requirements.txt" ]]; then
    echo "Installing Python packages (requirements.txt)..."
    pip3 install --upgrade pip
    pip3 install -r /workspaces/GAPMananas/requirements.txt || true
  fi
fi

# R packages (kept minimal) - run non-interactively
R -e "options(repos=list(CRAN='https://cloud.r-project.org')); if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager'); BiocManager::install(c('dada2','phyloseq'), ask=FALSE)" || true

echo "Post-create complete"
