#!/bin/bash
# Post-create script for GAPMananas bioinformatics pipeline
# Installs R packages, Python packages, and system dependencies

set -e

echo "=== Installing system dependencies ==="
sudo apt-get update
sudo apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    git

echo "=== Installing R packages ==="
# Install renv and BiocManager first
R -e "install.packages(c('renv', 'BiocManager'), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages for amplicon sequencing analysis
R -e "BiocManager::install(c('dada2', 'phyloseq', 'DECIPHER', 'decontam', 'phangorn'), ask=FALSE)"

# Additional packages for plotting and analysis
R -e "install.packages(c('ggplot2', 'tidyverse', 'vegan', 'ape', 'ampvis2'), repos='https://cloud.r-project.org/')"

echo "=== Installing Python packages ==="
pip3 install --upgrade pip
pip3 install -r /workspaces/GAPMananas/requirements.txt

# Install rpy2 inside the container where R is present (avoid CI build-time issues)
echo "=== Installing rpy2 (built against container R) ==="
# Force ABI mode during build to be more portable; limit cache to reduce image size
RPY2_CFFI_MODE=ABI pip3 install --no-cache-dir rpy2

echo "=== Restoring renv environment if available ==="
if [ -f "/workspaces/GAPMananas/renv.lock" ]; then
    cd /workspaces/GAPMananas
    R -e "renv::restore()"
fi

echo "=== Setup complete! ==="
echo "Jupyter Notebook: jupyter notebook --ip=0.0.0.0 --allow-root --no-browser"
echo "R Console: R"
