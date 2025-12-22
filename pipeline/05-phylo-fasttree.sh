#!/bin/bash

###############################################################################
# NucleiTaxa Pipeline Stage 05: Phylogenetic Tree Construction
# 
# FastTree 2 (ultra-fast maximum-likelihood tree inference)
# - 1000x faster than UPGMA or neighbor-joining
# - Better accuracy for large alignments (10K+ sequences)
# - Suitable for microbiome ASV trees
#
# Input:  Representative ASV sequences from Stage 03
# Output: Phylogenetic tree (Newick format), rooted tree
###############################################################################

set -euo pipefail

# Default parameters
INPUT_DIR=""
OUTPUT_DIR=""
PROFILE="default"
JOBS=4
CONFIG_FILE=""
ANALYTICS_PORT=8888
LOG_FILE=""

# FastTree parameters
SEQUENCE_TYPE="DNA"  # DNA, AA (amino acid)
OUTGROUP=""  # Optional outgroup sequence ID for rooting
BOOT_REPLICATES=0  # 0 = no bootstrap, >0 = number of replicates

# Logging
log_info() {
    echo "[INFO] $*" | tee -a "$LOG_FILE"
}

log_success() {
    echo "[✓] $*" | tee -a "$LOG_FILE"
}

log_error() {
    echo "[ERROR] $*" | tee -a "$LOG_FILE"
}

log_warning() {
    echo "[⚠] $*" | tee -a "$LOG_FILE"
}

emit_metric() {
    local key="$1"
    local value="$2"
    
    if command -v nc &>/dev/null; then
        echo "$key=$value" | nc -w1 localhost "$ANALYTICS_PORT" 2>/dev/null || true
    fi
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)       INPUT_DIR="$2"; shift 2 ;;
        --output-dir)      OUTPUT_DIR="$2"; shift 2 ;;
        --profile)         PROFILE="$2"; shift 2 ;;
        --jobs)            JOBS="$2"; shift 2 ;;
        --config)          CONFIG_FILE="$2"; shift 2 ;;
        --analytics-port)  ANALYTICS_PORT="$2"; shift 2 ;;
        *)                 log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Usage: $0 --input-dir <path> --output-dir <path>"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/05-phylo"
LOG_FILE="$OUTPUT_DIR/05-phylo/phylo.log"

log_info "=== NucleiTaxa Stage 05: Phylogenetic Tree (FastTree 2) ==="
log_info "Input:  $OUTPUT_DIR/03-chimera"
log_info "Output: $OUTPUT_DIR/05-phylo"
log_info "Profile: $PROFILE"

# Load profile settings
if [[ -f "$CONFIG_FILE" ]]; then
    log_info "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE" || true
fi

# Verify FastTree is installed
if ! command -v fasttree &>/dev/null && ! command -v FastTree &>/dev/null; then
    log_error "FastTree 2 not found. Install with: apt install fasttree"
    exit 1
fi

# Get FastTree version
FT_CMD=$(command -v fasttree || command -v FastTree)
FASTTREE_VERSION=$("$FT_CMD" 2>&1 | grep -oP 'FastTree.*' | head -1 || echo "FastTree 2")
log_success "$FASTTREE_VERSION verified"

# Check input file
ASV_FASTA="$OUTPUT_DIR/03-chimera/nonchimeras_denovo.fasta"
if [[ ! -f "$ASV_FASTA" ]]; then
    log_error "ASV sequences not found: $ASV_FASTA"
    exit 1
fi

SEQ_COUNT=$(grep -c "^>" "$ASV_FASTA")
log_success "Sequences to align and tree: $SEQ_COUNT"

# Stage 1: Multiple Sequence Alignment (if needed)
log_info "=== Sequence Alignment (if needed) ==="

# For moderate-sized datasets (10K-100K seqs), MUSCLE/MAFFT is too slow
# FastTree can use unaligned sequences with -nt (nucleotide) flag
# Alternative: use alignment acceleration (clustalo)

if command -v clustalo &>/dev/null && [[ $SEQ_COUNT -lt 50000 ]]; then
    log_info "Using Clustal Omega for alignment..."
    
    clustalo \
        -i "$ASV_FASTA" \
        -o "$OUTPUT_DIR/05-phylo/asv_aligned.fasta" \
        -t DNA \
        --threads="$JOBS" \
        --output-order=input-order \
        2>&1 | tee -a "$LOG_FILE"
    
    ALIGNMENT_INPUT="$OUTPUT_DIR/05-phylo/asv_aligned.fasta"
    log_success "Alignment complete"
else
    log_info "Using unaligned input (FastTree -nt mode)"
    ALIGNMENT_INPUT="$ASV_FASTA"
fi

# Stage 2: Build phylogenetic tree with FastTree
log_info "=== Building Maximum-Likelihood Tree ==="

# FastTree flags:
# -nt:           Treat input as DNA (unaligned)
# -n:            Maximum likelihood tree (default)
# -gtr:          Use GTR model (more accurate than JC69)
# -gamma:        Model gamma distribution of rates (accounting for site-to-site variation)
# -quiet:        Minimal output
# Additional performance: -mlnni 4 (ML nearest neighbor interchanges)

FASTTREE_OPTS="-nt -gtr -gamma -quiet"

log_info "FastTree options: $FASTTREE_OPTS"

"$FT_CMD" $FASTTREE_OPTS \
    "$ALIGNMENT_INPUT" \
    > "$OUTPUT_DIR/05-phylo/asv_tree_unrooted.nwk" \
    2>&1 | tee -a "$LOG_FILE"

if [[ ! -f "$OUTPUT_DIR/05-phylo/asv_tree_unrooted.nwk" ]]; then
    log_error "FastTree failed to generate tree"
    exit 1
fi

log_success "Unrooted tree generated"

# Stage 3: Root the tree (optional, midpoint rooting as default)
log_info "=== Rooting Phylogenetic Tree ==="

# Use R's ape package to root tree at midpoint
# Midpoint rooting is unbiased (works when outgroup unknown)

cat > "$OUTPUT_DIR/05-phylo/root_tree.R" << 'EOFR'
#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly = TRUE)
unrooted_tree_file <- args[1]
output_dir <- args[2]
outgroup <- if (length(args) >= 3 && args[3] != "") args[3] else NULL

cat("Loading tree from:", unrooted_tree_file, "\n")
tree <- read.tree(unrooted_tree_file)

cat("Unrooted tree: ", length(tree$tip.label), "tips\n")

# Root tree at midpoint
tree_rooted <- midpoint.root(tree)

cat("Rooting complete. Rooted tree has", length(tree_rooted$tip.label), "tips\n")

# Save rooted tree
output_file <- file.path(output_dir, "asv_tree_rooted.nwk")
write.tree(tree_rooted, output_file)
cat("Rooted tree saved to:", output_file, "\n")

# Generate tree statistics
cat("\n=== Tree Statistics ===\n")
cat("Number of tips:", length(tree_rooted$tip.label), "\n")
cat("Number of internal nodes:", tree_rooted$Nnode, "\n")
cat("Tree height (max distance from root):", max(dist.nodes(tree_rooted)[1, ]), "\n")

# Save statistics
stats_file <- file.path(output_dir, "tree_stats.txt")
cat("Tree Statistics\n", file = stats_file)
cat("Number of sequences:", length(tree_rooted$tip.label), "\n", file = stats_file, append = TRUE)
cat("Number of internal nodes:", tree_rooted$Nnode, "\n", file = stats_file, append = TRUE)
cat("Tree height:", max(dist.nodes(tree_rooted)[1, ]), "\n", file = stats_file, append = TRUE)

cat("Tree statistics saved\n")
EOFR

if command -v Rscript &>/dev/null; then
    Rscript "$OUTPUT_DIR/05-phylo/root_tree.R" \
        "$OUTPUT_DIR/05-phylo/asv_tree_unrooted.nwk" \
        "$OUTPUT_DIR/05-phylo" \
        "$OUTGROUP" \
        2>&1 | tee -a "$LOG_FILE" || log_warning "Tree rooting failed (using unrooted)"
    
    if [[ -f "$OUTPUT_DIR/05-phylo/asv_tree_rooted.nwk" ]]; then
        FINAL_TREE="$OUTPUT_DIR/05-phylo/asv_tree_rooted.nwk"
    else
        FINAL_TREE="$OUTPUT_DIR/05-phylo/asv_tree_unrooted.nwk"
    fi
else
    log_warning "R not available, cannot root tree. Using unrooted tree."
    FINAL_TREE="$OUTPUT_DIR/05-phylo/asv_tree_unrooted.nwk"
fi

log_success "Final tree: $FINAL_TREE"

# Stage 4: Validate tree structure
log_info "=== Validating Tree Structure ==="

# Check Newick format is valid
if grep -q "^[a-zA-Z]" "$FINAL_TREE"; then
    LEAF_COUNT=$(grep -o "[a-zA-Z0-9_]*:" "$FINAL_TREE" | wc -l)
    log_success "Tree validation: $LEAF_COUNT leaves detected"
else
    log_error "Tree format may be invalid"
fi

# Stage 5: Generate phylogenetic diversity metrics
log_info "=== Computing Phylogenetic Diversity Metrics ==="

cat > "$OUTPUT_DIR/05-phylo/phylo_metrics.R" << 'EOFR'
#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
seqtab_file <- args[2]
output_dir <- args[3]

cat("Loading phylogenetic tree:", tree_file, "\n")
tree <- read.tree(tree_file)

cat("Loading sequence table:", seqtab_file, "\n")
seqtab <- read.table(seqtab_file, sep = "\t", header = TRUE, row.names = 1)

cat("Computing phylogenetic metrics...\n")

# Phylogenetic diversity = sum of branch lengths
tree_length <- sum(tree$edge.length)
cat("Total tree length:", round(tree_length, 2), "\n")

# Average distance between tips
distances <- dist.nodes(tree)
tip_distances <- distances[1:length(tree$tip.label), 1:length(tree$tip.label)]
avg_dist <- mean(tip_distances[upper.tri(tip_distances)])
cat("Average pairwise distance:", round(avg_dist, 4), "\n")

# Save metrics
metrics_file <- file.path(output_dir, "phylo_metrics.txt")
cat("Phylogenetic Metrics\n", file = metrics_file)
cat("Total tree length:", tree_length, "\n", file = metrics_file, append = TRUE)
cat("Average pairwise distance:", avg_dist, "\n", file = metrics_file, append = TRUE)
cat("Metrics saved\n")
EOFR

if command -v Rscript &>/dev/null; then
    Rscript "$OUTPUT_DIR/05-phylo/phylo_metrics.R" \
        "$FINAL_TREE" \
        "$OUTPUT_DIR/03-chimera/seqtab_nochim.txt" \
        "$OUTPUT_DIR/05-phylo" \
        2>&1 | tee -a "$LOG_FILE" || log_warning "Metrics computation failed (non-critical)"
fi

# Generate comprehensive report
cat > "$OUTPUT_DIR/05-phylo/PHYLO_REPORT.txt" << EOF
=== Phylogenetic Analysis Report ===
Timestamp: $(date)
Profile: $PROFILE

METHOD: FastTree 2 (Maximum Likelihood)
Algorithm: GTR model with gamma-distributed rate variation
Sequences analyzed: $SEQ_COUNT ASVs
Alignment: $([ -f "$OUTPUT_DIR/05-phylo/asv_aligned.fasta" ] && echo "Clustal Omega" || echo "FastTree -nt (unaligned)")

TREE PROPERTIES:
  Tree file: $(basename "$FINAL_TREE")
  Rooted: $([ "$FINAL_TREE" = "$OUTPUT_DIR/05-phylo/asv_tree_rooted.nwk" ] && echo "Yes (midpoint)" || echo "No (unrooted)")

OUTPUTS:
  - $(basename "$FINAL_TREE") (primary phylogenetic tree)
  - tree_stats.txt (tree statistics)
  - phylo_metrics.txt (phylogenetic diversity)
  
NEXT STAGE:
  Proceed to 06-krona-viz.sh for interactive taxonomy visualization

INTERPRETATION:
  Tree branch lengths represent evolutionary distance
  Use for community structure analysis, phylogenetic diversity metrics
  Compatible with QIIME2 / PhyloSeq for downstream analysis

REFERENCES:
  - FastTree: Price et al., 2010 (PLoS ONE)
  - GTR model: Tavare, 1986
EOF

log_success "Phylogenetic report written"

# Emit metrics
emit_metric "tree_generated" "1"
emit_metric "sequence_count_in_tree" "$SEQ_COUNT"

log_success "=== Stage 05 Complete ==="
