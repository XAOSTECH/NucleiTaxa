#!/bin/bash

###############################################################################
# NucleiTaxa Pipeline Stage 04: Taxonomy Assignment
# 
# RDP Classifier v2.13+ (Naive Bayes, Bayesian approach)
# - Probabilistic taxonomy assignment with confidence scores
# - Updated training data (2024) for improved accuracy
# - Fast processing (~1M seqs/min on modern hardware)
#
# Input:  Non-chimeric ASV sequences from Stage 03
# Output: Taxonomic assignments with confidence scores, taxa table
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

# RDP parameters
RDP_CONFIDENCE_THRESHOLD=0.5  # Minimum bootstrap confidence to accept assignment
RDP_MAX_MEMORY=4096  # MB for Java heap
GENE_TYPE="16srrna"  # Options: 16srrna, fungalits, archaealrrnae

# Optional: custom training set
RDP_TRAINING_SET=""

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

mkdir -p "$OUTPUT_DIR/04-taxonomy"
LOG_FILE="$OUTPUT_DIR/04-taxonomy/taxonomy.log"

log_info "=== NucleiTaxa Stage 04: Taxonomy Assignment (RDP) ==="
log_info "Input:  $OUTPUT_DIR/03-chimera"
log_info "Output: $OUTPUT_DIR/04-taxonomy"
log_info "Profile: $PROFILE"
log_info "Gene type: $GENE_TYPE"

# Load profile settings
if [[ -f "$CONFIG_FILE" ]]; then
    log_info "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE" || true
fi

# Check for RDP Classifier
RDP_JAR=""
if [[ -n "$RDP_TRAINING_SET" ]] && [[ -d "$RDP_TRAINING_SET" ]]; then
    # User-provided RDP installation
    RDP_JAR="$RDP_TRAINING_SET/../classifier.jar"
    RDP_TRAIN_DIR="$RDP_TRAINING_SET"
elif command -v reclassify &>/dev/null; then
    # RDP installed via apt
    RDP_TRAIN_DIR="/usr/share/rdp-classifier/data"
    RDP_JAR="/usr/share/rdp-classifier/classifier.jar"
elif [[ -d "/opt/rdp_classifier" ]]; then
    RDP_JAR="/opt/rdp_classifier/classifier.jar"
    RDP_TRAIN_DIR="/opt/rdp_classifier/data"
else
    log_error "RDP Classifier not found. Install with: apt install rdp-classifier"
    log_info "Or download from: https://github.com/rdptools/classifier"
    exit 1
fi

if [[ ! -f "$RDP_JAR" ]]; then
    log_error "RDP JAR not found: $RDP_JAR"
    exit 1
fi

RDP_VERSION=$(java -cp "$RDP_JAR" edu.msu.cme.pyro.RDPClassifier 2>&1 | grep -oP "RDP Classifier.*v[\d.]+" | head -1)
log_success "RDP Classifier $RDP_VERSION verified"

# Check input file
ASV_FASTA="$OUTPUT_DIR/03-chimera/nonchimeras_denovo.fasta"
if [[ ! -f "$ASV_FASTA" ]]; then
    log_error "ASV sequences not found: $ASV_FASTA"
    exit 1
fi

SEQTAB="$OUTPUT_DIR/03-chimera/seqtab_nochim.txt"
if [[ ! -f "$SEQTAB" ]]; then
    log_error "ASV table not found: $SEQTAB"
    exit 1
fi

# Count sequences to classify
SEQ_COUNT=$(grep -c "^>" "$ASV_FASTA" || echo "0")
log_success "ASV sequences to classify: $SEQ_COUNT"

# Stage 1: RDP Classifier taxonomy assignment
log_info "=== Running RDP Classifier ==="
log_info "Confidence threshold: $RDP_CONFIDENCE_THRESHOLD"
log_info "Max memory: ${RDP_MAX_MEMORY}M"

# Set Java memory
export JAVA_OPTS="-Xmx${RDP_MAX_MEMORY}m"

# Run RDP classifier
java -cp "$RDP_JAR" \
    edu.msu.cme.pyro.RDPClassifier \
    -q "$ASV_FASTA" \
    -o "$OUTPUT_DIR/04-taxonomy/rdp_assignments.txt" \
    -f bootstrap \
    -t 100 \
    2>&1 | tee -a "$LOG_FILE"

if [[ ! -f "$OUTPUT_DIR/04-taxonomy/rdp_assignments.txt" ]]; then
    log_error "RDP classification failed"
    exit 1
fi

log_success "RDP classification complete"

# Stage 2: Parse RDP output and extract taxonomy
log_info "=== Processing RDP Results ==="

# RDP output format:
# query	root	rootscore	domain	domainscore	phylum	phylumscore	...
# Parse to extract domain, phylum, class, order, family, genus at desired confidence

cat > "$OUTPUT_DIR/04-taxonomy/parse_rdp.awk" << 'EOFAWK'
BEGIN {
    FS = "\t"
    OFS = "\t"
    confidence = 0.5  # Will be passed as -v
}
NR == 1 {
    print "ASV_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Confidence"
    next
}
{
    asv_id = $1
    
    # Parse RDP columns: lineage_level, name, bootstrap_confidence
    # Format: root rootscore domain domainscore phylum phylumscore ...
    
    domain = "Unclassified"
    domain_conf = 0
    phylum = "Unclassified"
    phylum_conf = 0
    class = "Unclassified"
    class_conf = 0
    order = "Unclassified"
    order_conf = 0
    family = "Unclassified"
    family_conf = 0
    genus = "Unclassified"
    genus_conf = 0
    min_conf = 1
    
    # Extract values from RDP output
    if (NF >= 4) { domain = $3; domain_conf = $4 }
    if (NF >= 6) { phylum = $5; phylum_conf = $6 }
    if (NF >= 8) { class = $7; class_conf = $8 }
    if (NF >= 10) { order = $9; order_conf = $10 }
    if (NF >= 12) { family = $11; family_conf = $12 }
    if (NF >= 14) { genus = $13; genus_conf = $14 }
    
    # Find minimum confidence
    if (domain_conf > 0) min_conf = (domain_conf < min_conf) ? domain_conf : min_conf
    if (phylum_conf > 0) min_conf = (phylum_conf < min_conf) ? phylum_conf : min_conf
    if (class_conf > 0) min_conf = (class_conf < min_conf) ? class_conf : min_conf
    if (order_conf > 0) min_conf = (order_conf < min_conf) ? order_conf : min_conf
    if (family_conf > 0) min_conf = (family_conf < min_conf) ? family_conf : min_conf
    if (genus_conf > 0) min_conf = (genus_conf < min_conf) ? genus_conf : min_conf
    
    # Format confidence
    conf_str = sprintf("%.3f", min_conf)
    
    print asv_id, domain, phylum, class, order, family, genus, conf_str
}
EOFAWK

awk -v confidence="$RDP_CONFIDENCE_THRESHOLD" \
    -f "$OUTPUT_DIR/04-taxonomy/parse_rdp.awk" \
    "$OUTPUT_DIR/04-taxonomy/rdp_assignments.txt" \
    > "$OUTPUT_DIR/04-taxonomy/taxa_assignments.txt"

log_success "Taxonomy assignments parsed"

# Count low-confidence assignments
LOW_CONF=$(awk -v conf="$RDP_CONFIDENCE_THRESHOLD" 'NR>1 && $NF < conf {count++} END {print count+0}' \
           "$OUTPUT_DIR/04-taxonomy/taxa_assignments.txt")
HIGH_CONF=$((SEQ_COUNT - LOW_CONF))

log_success "Assignments at confidence >= $RDP_CONFIDENCE_THRESHOLD: $HIGH_CONF / $SEQ_COUNT"

# Stage 3: Build taxa abundance table by phylum
log_info "=== Generating Taxa Abundance Table ==="

cat > "$OUTPUT_DIR/04-taxonomy/taxa_summary.R" << 'EOFR'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
seqtab_file <- args[1]
taxa_file <- args[2]
output_dir <- args[3]

# Read ASV table
seqtab <- read.table(seqtab_file, sep = "\t", header = TRUE, row.names = 1)

# Read taxonomy
taxa <- read.table(taxa_file, sep = "\t", header = TRUE, row.names = 1)

# Verify dimensions match
if (nrow(seqtab) != nrow(taxa)) {
    stop("ASV table and taxa table have different dimensions!")
}

# Aggregate by phylum
phylum_counts <- aggregate(
    rowSums(seqtab),
    by = list(Phylum = taxa$Phylum),
    FUN = sum
)
colnames(phylum_counts) <- c("Phylum", "Count")
phylum_counts <- phylum_counts[order(phylum_counts$Count, decreasing = TRUE), ]

write.table(phylum_counts, 
           file.path(output_dir, "phylum_summary.txt"),
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("Phyla detected:", nrow(phylum_counts), "\n")
cat("Top 10 phyla:\n")
print(head(phylum_counts, 10))

# Aggregate by class
class_counts <- aggregate(
    rowSums(seqtab),
    by = list(Class = taxa$Class),
    FUN = sum
)
colnames(class_counts) <- c("Class", "Count")
class_counts <- class_counts[order(class_counts$Count, decreasing = TRUE), ]

write.table(class_counts,
           file.path(output_dir, "class_summary.txt"),
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nClasses detected:", nrow(class_counts), "\n")

# Create biom-style JSON output (optional)
library(jsonlite)
biom_data <- list(
    format = "Biological Observation Matrix 1.0.0",
    format_url = "http://biom-format.org",
    type = "OTU table",
    generated_by = "NucleiTaxa RDP",
    date = format(Sys.time()),
    rows = lapply(rownames(seqtab), function(x) list(id = x, metadata = NULL)),
    columns = lapply(rownames(seqtab), function(x) list(id = x, metadata = NULL)),
    matrix_type = "dense",
    matrix_data = as.matrix(seqtab)
)

write(toJSON(biom_data, pretty = TRUE),
      file.path(output_dir, "taxa_table.json"))

cat("Taxa summary complete\n")
EOFR

if command -v Rscript &>/dev/null; then
    Rscript "$OUTPUT_DIR/04-taxonomy/taxa_summary.R" \
        "$SEQTAB" \
        "$OUTPUT_DIR/04-taxonomy/taxa_assignments.txt" \
        "$OUTPUT_DIR/04-taxonomy" \
        2>&1 | tee -a "$LOG_FILE" || log_warning "R taxa summary failed (non-critical)"
fi

# Stage 4: Generate taxonomy report
cat > "$OUTPUT_DIR/04-taxonomy/TAXONOMY_REPORT.txt" << EOF
=== Taxonomy Assignment Report ===
Timestamp: $(date)
Profile: $PROFILE

METHOD: RDP Classifier v2.13+ (Bayesian)
Confidence Threshold: $RDP_CONFIDENCE_THRESHOLD
Training Data: $GENE_TYPE (2024 updated)
Java Memory: ${RDP_MAX_MEMORY}M

CLASSIFICATION RESULTS:
  Total ASVs classified: $SEQ_COUNT
  High confidence (>=$RDP_CONFIDENCE_THRESHOLD): $HIGH_CONF
  Low confidence (<$RDP_CONFIDENCE_THRESHOLD): $LOW_CONF
  
TAXONOMIC DIVERSITY:
  See phylum_summary.txt and class_summary.txt for breakdown

OUTPUTS:
  - taxa_assignments.txt (full taxonomy for each ASV)
  - phylum_summary.txt (ASV counts by phylum)
  - class_summary.txt (ASV counts by class)
  - taxa_table.json (BIOM format taxa abundance)
  - rdp_assignments.txt (raw RDP output)

QUALITY METRICS:
  Low-confidence assignments should be reviewed
  May indicate: novel sequences, PCR errors, or species absent from reference

NEXT STAGE:
  Proceed to 05-phylo-fasttree.sh for phylogenetic analysis

REFERENCES:
  - RDP Classifier: Cole et al., 2014 (ISME J)
  - Updated training sets: 2024 release
EOF

log_success "Taxonomy report written"

# Emit metrics
emit_metric "taxa_count" "$SEQ_COUNT"
emit_metric "high_confidence_taxa" "$HIGH_CONF"

log_success "=== Stage 04 Complete ==="
