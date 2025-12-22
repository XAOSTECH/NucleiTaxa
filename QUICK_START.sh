#!/bin/bash

###############################################################################
# NucleiTaxa QUICK START REFERENCE
# 
# Copy-paste recipes for common analyses
###############################################################################

# ============================================================================
# BASIC EXECUTION
# ============================================================================

# Minimal command (everything automated)
./bin/nucleitaxa --forward input_R1.fastq.gz --reverse input_R2.fastq.gz

# With custom output directory
./bin/nucleitaxa \
    --forward input_R1.fastq.gz \
    --reverse input_R2.fastq.gz \
    --output my_results

# ============================================================================
# WITH REAL-TIME ANALYTICS
# ============================================================================

# Terminal 1: Start analytics server
cd /path/to/NucleiTaxa
./analytics/server/nucleitaxa-server &

# Terminal 2: Run pipeline with analytics enabled
./bin/nucleitaxa \
    --forward input_R1.fastq.gz \
    --reverse input_R2.fastq.gz \
    --analytics-live localhost:8888

# Then open browser: http://localhost:8888

# ============================================================================
# DATASET-SPECIFIC PROFILES
# ============================================================================

# 16S rRNA (default, bacterial/archaeal)
./bin/nucleitaxa \
    --profile 16s \
    --forward sample_R1.fastq.gz \
    --reverse sample_R2.fastq.gz

# ITS (fungal Internal Transcribed Spacer)
./bin/nucleitaxa \
    --profile its \
    --forward sample_R1.fastq.gz \
    --reverse sample_R2.fastq.gz

# ============================================================================
# PARALLEL PROCESSING & MEMORY
# ============================================================================

# Use 16 CPU cores (default: 4)
./bin/nucleitaxa \
    --forward R1.fastq.gz \
    --reverse R2.fastq.gz \
    --jobs 16

# For large datasets (increase RDP memory)
export RDP_MAX_MEMORY=8192  # 8GB
./bin/nucleitaxa --forward R1.fastq.gz --reverse R2.fastq.gz

# ============================================================================
# CUSTOM CONFIGURATION
# ============================================================================

# Create configuration file
cat > my_settings.cfg << 'EOF'
TRUNCATE_LENGTH_F=250
TRUNCATE_LENGTH_R=200
MIN_QUALITY=20
MAX_EE_F=2
MAX_EE_R=2
MIN_OVERLAP=12
GENE_TYPE="16srrna"
UCHIME_THRESHOLD=0.85
RDP_CONFIDENCE_THRESHOLD=0.5
EOF

# Use configuration file
./bin/nucleitaxa \
    --forward R1.fastq.gz \
    --reverse R2.fastq.gz \
    --config my_settings.cfg

# ============================================================================
# TROUBLESHOOTING & VALIDATION
# ============================================================================

# Dry-run (show commands without executing)
./bin/nucleitaxa \
    --forward R1.fastq.gz \
    --reverse R2.fastq.gz \
    --dry-run

# Resume from stage 03 (if pipeline interrupted)
./bin/nucleitaxa \
    --resume-from 03 \
    --output results

# Run only stage 04 (taxonomy)
bash pipeline/04-taxonomy-rdp.sh \
    --input-dir results \
    --output-dir results

# Validate pipeline structure (with mock data)
bash test-suite.sh

# ============================================================================
# VIEWING RESULTS
# ============================================================================

# Open interactive taxonomy visualization
# (Replace with your file path)
open results/06-viz/taxa_krona.html

# View ASV abundance table
cat results/03-chimera/seqtab_nochim.txt | head

# Check taxonomy assignments
head results/04-taxonomy/taxa_assignments.txt

# View phylogenetic tree
cat results/05-phylo/asv_tree_rooted.nwk

# ============================================================================
# BATCH PROCESSING MULTIPLE SAMPLES
# ============================================================================

# Process multiple sample pairs in parallel
for sample in sample1 sample2 sample3; do
    ./bin/nucleitaxa \
        --forward data/${sample}_R1.fastq.gz \
        --reverse data/${sample}_R2.fastq.gz \
        --output results/${sample} \
        --jobs 4 &
done
wait  # Wait for all to complete

# ============================================================================
# DOWNSTREAM ANALYSIS INTEGRATION
# ============================================================================

# Export to QIIME2 format
# (After pipeline completion)
# Copy results to QIIME2 import directory:
cp results/03-chimera/seqtab_nochim.txt qiime2_import/
cp results/04-taxonomy/taxa_assignments.txt qiime2_import/
cp results/05-phylo/asv_tree_rooted.nwk qiime2_import/

# Import into PhyloSeq (R)
# R code:
# library(phyloseq)
# seqtab <- read.table("results/03-chimera/seqtab_nochim.txt", sep="\t", header=T, row.names=1)
# tax <- read.table("results/04-taxonomy/taxa_assignments.txt", sep="\t", header=T, row.names=1)
# tree <- read_tree("results/05-phylo/asv_tree_rooted.nwk")
# ps <- phyloseq(otu_table(seqtab, taxa_are_rows=F), tax_table(as.matrix(tax)), tree)

# ============================================================================
# INSTALLATION QUICK CHECK
# ============================================================================

# Verify all required tools
echo "Checking dependencies..."

command -v R > /dev/null && echo "✓ R" || echo "✗ R (install: apt install r-base)"
Rscript -e "library(dada2)" 2>/dev/null && echo "✓ DADA2" || echo "✗ DADA2"
command -v vsearch > /dev/null && echo "✓ VSEARCH" || echo "✗ VSEARCH"
command -v fasttree > /dev/null && echo "✓ FastTree" || echo "✗ FastTree"
command -v java > /dev/null && echo "✓ Java (for RDP)" || echo "✗ Java"
command -v ktImportTaxonomy > /dev/null && echo "✓ Krona" || echo "✗ Krona"

# ============================================================================
# PERFORMANCE TUNING
# ============================================================================

# For slow systems (limit memory)
export RDP_MAX_MEMORY=2048  # 2GB
./bin/nucleitaxa --jobs 2 --forward R1.fastq.gz --reverse R2.fastq.gz

# For fast systems (maximize parallelism)
export RDP_MAX_MEMORY=16384  # 16GB
./bin/nucleitaxa --jobs 32 --forward R1.fastq.gz --reverse R2.fastq.gz

# ============================================================================
# DEBUGGING
# ============================================================================

# Enable verbose logging (check nucleitaxa script for log redirection)
set -x  # Show all commands
./bin/nucleitaxa --forward R1.fastq.gz --reverse R2.fastq.gz
set +x  # Turn off verbose

# Check individual stage logs
tail -f results/01-preprocess/preprocess.log
tail -f results/02-denoise/denoise.log
tail -f results/03-chimera/chimera.log
tail -f results/04-taxonomy/taxonomy.log
tail -f results/05-phylo/phylo.log
tail -f results/06-viz/visualization.log

# ============================================================================
# COMMON ERRORS & SOLUTIONS
# ============================================================================

# Error: "DADA2 R package not found"
# Solution:
Rscript -e "install.packages('BiocManager'); BiocManager::install('dada2')"

# Error: "VSEARCH not found"
# Solution:
sudo apt install -y vsearch

# Error: "RDP Classifier failed"
# Solution:
export JAVA_OPTS="-Xmx4096m"  # Increase Java heap
./bin/nucleitaxa --forward R1.fastq.gz --reverse R2.fastq.gz

# Error: "No R1 filtered FASTQ files found"
# Solution:
# Check that 01-preprocess stage ran successfully
ls -la results/01-preprocess/
# Verify input files exist:
ls -la input_R1.fastq.gz input_R2.fastq.gz

# ============================================================================
# SUMMARY OF KEY FILES
# ============================================================================

# Main orchestrator
# ./bin/nucleitaxa

# Pipeline stages
# ./pipeline/01-preprocess.sh          → BBTools QC
# ./pipeline/02-denoise-dada2.sh       → DADA2 ASV inference
# ./pipeline/03-chimera-vsearch.sh     → VSEARCH UCHIME
# ./pipeline/04-taxonomy-rdp.sh        → RDP Classifier
# ./pipeline/05-phylo-fasttree.sh      → FastTree phylogenetics
# ./pipeline/06-krona-viz.sh           → Krona visualization

# Analytics
# ./analytics/server/nucleitaxa-server.cpp   → WebSocket backend
# ./analytics/web/app.js                     → Frontend client
# ./analytics/web/index.html                 → Dashboard UI
# ./analytics/web/styles.css                 → Styling

# Test & documentation
# ./test-suite.sh                     → Validation
# ./README_BASH_REFACTOR_UPDATED.md   → Full documentation

echo "Quick reference complete. For detailed info, see README_BASH_REFACTOR_UPDATED.md"
