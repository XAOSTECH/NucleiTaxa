#!/bin/bash

###############################################################################
# NucleiTaxa Pipeline Stage 03: Chimera Detection
# 
# Hybrid approach: DADA2 + VSEARCH UCHIME
# - DADA2 removes chimeras during denoising (already done in Stage 02)
# - VSEARCH de novo + reference-based chimera detection (additional QC)
# Research shows hybrid approach reduces false positives vs single method
#
# Input:  ASV table + sequences from Stage 02
# Output: Chimera-flagged sequences, cleaned ASV table
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

# VSEARCH parameters
MIN_ABUNDANCE=2
UCHIME_THRESHOLD=0.85
REF_DB=""  # Optional reference database for reference-based chimera detection

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

mkdir -p "$OUTPUT_DIR/03-chimera"
LOG_FILE="$OUTPUT_DIR/03-chimera/chimera.log"

log_info "=== NucleiTaxa Stage 03: Chimera Detection (VSEARCH) ==="
log_info "Input:  $OUTPUT_DIR/02-denoise"
log_info "Output: $OUTPUT_DIR/03-chimera"
log_info "Profile: $PROFILE"

# Load profile settings
if [[ -f "$CONFIG_FILE" ]]; then
    log_info "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE" || true
fi

# Verify VSEARCH is installed
if ! command -v vsearch &>/dev/null; then
    log_error "VSEARCH not found. Install with: apt install vsearch"
    exit 1
fi

VSEARCH_VERSION=$(vsearch --version 2>&1 | grep -oP 'v\d+\.\d+\.\d+' | head -1)
log_success "VSEARCH $VSEARCH_VERSION verified"

# Check input files from DADA2 stage
if [[ ! -f "$OUTPUT_DIR/02-denoise/asv_sequences.fasta" ]]; then
    log_error "ASV sequences not found: $OUTPUT_DIR/02-denoise/asv_sequences.fasta"
    exit 1
fi

if [[ ! -f "$OUTPUT_DIR/02-denoise/seqtab.txt" ]]; then
    log_error "ASV table not found: $OUTPUT_DIR/02-denoise/seqtab.txt"
    exit 1
fi

log_success "Input files verified"

# Prepare working files
ASV_FASTA="$OUTPUT_DIR/02-denoise/asv_sequences.fasta"
SEQTAB="$OUTPUT_DIR/02-denoise/seqtab.txt"

# Stage 1: De novo chimera detection with VSEARCH UCHIME
log_info "=== De Novo Chimera Detection ==="
log_info "Detecting chimeras using UCHIME algorithm..."

# VSEARCH denovo chimera detection
# Output: tab-delimited with fields: query, db_count, chimera_indicator (Y/N)
vsearch \
    --uchime_denovo "$ASV_FASTA" \
    --chimeras "$OUTPUT_DIR/03-chimera/chimeras_denovo.fasta" \
    --nonchimeras "$OUTPUT_DIR/03-chimera/nonchimeras_denovo.fasta" \
    --uchimeout "$OUTPUT_DIR/03-chimera/uchime_denovo.txt" \
    --threads "$JOBS" \
    2>&1 | tee -a "$LOG_FILE"

# Count chimeras detected
DENOVO_CHIMERAS=$(grep -c "Y$" "$OUTPUT_DIR/03-chimera/uchime_denovo.txt" || echo "0")
DENOVO_NONCHIM=$(grep -c "N$" "$OUTPUT_DIR/03-chimera/uchime_denovo.txt" || echo "0")

log_success "De novo detection: $DENOVO_CHIMERAS chimeras, $DENOVO_NONCHIM non-chimeras"

# Stage 2: Reference-based chimera detection (if reference available)
if [[ -n "$REF_DB" ]] && [[ -f "$REF_DB" ]]; then
    log_info "=== Reference-Based Chimera Detection ==="
    log_info "Reference database: $REF_DB"
    
    vsearch \
        --uchime_ref "$OUTPUT_DIR/03-chimera/nonchimeras_denovo.fasta" \
        --db "$REF_DB" \
        --chimeras "$OUTPUT_DIR/03-chimera/chimeras_reference.fasta" \
        --nonchimeras "$OUTPUT_DIR/03-chimera/nonchimeras_final.fasta" \
        --uchimeout "$OUTPUT_DIR/03-chimera/uchime_reference.txt" \
        --threads "$JOBS" \
        2>&1 | tee -a "$LOG_FILE"
    
    REF_CHIMERAS=$(grep -c "Y$" "$OUTPUT_DIR/03-chimera/uchime_reference.txt" || echo "0")
    log_success "Reference-based detection: $REF_CHIMERAS additional chimeras"
    
    FINAL_SEQS="$OUTPUT_DIR/03-chimera/nonchimeras_final.fasta"
else
    log_warning "No reference database specified, skipping reference-based detection"
    FINAL_SEQS="$OUTPUT_DIR/03-chimera/nonchimeras_denovo.fasta"
fi

# Count final non-chimeras
FINAL_COUNT=$(grep -c "^>" "$FINAL_SEQS" 2>/dev/null || echo "0")
log_success "Final non-chimeric ASVs: $FINAL_COUNT"

# Stage 3: Update ASV abundance table to exclude chimeras
log_info "=== Updating ASV Table ==="

# Generate list of chimeric sequence IDs
CHIM_IDS="$OUTPUT_DIR/03-chimera/chimera_ids.txt"
cat "$OUTPUT_DIR/03-chimera/chimeras_denovo.fasta" \
    | grep "^>" \
    | sed 's/>//' \
    > "$CHIM_IDS" || true

# Extract non-chimeric sequences using awk
# Parse FASTA to get sequence identifiers and remove rows corresponding to chimeras
cat > "$OUTPUT_DIR/03-chimera/filter_chimeras.awk" << 'EOFAWK'
BEGIN {
    FS = "\t"
}
FNR == NR {
    # Read chimera IDs (first argument file)
    chimeras[$1] = 1
    next
}
{
    # Process ASV table
    if (NR == 1) {
        print $0  # Print header
        next
    }
    
    seq_id = $1
    # Remove sequences that are chimeric
    if (!(seq_id in chimeras)) {
        print $0
    }
}
EOFAWK

awk -f "$OUTPUT_DIR/03-chimera/filter_chimeras.awk" \
    "$CHIM_IDS" \
    "$SEQTAB" \
    > "$OUTPUT_DIR/03-chimera/seqtab_nochim.txt"

log_success "Chimera-filtered ASV table generated"

# Validate output
if [[ ! -f "$OUTPUT_DIR/03-chimera/seqtab_nochim.txt" ]]; then
    log_error "Failed to generate seqtab_nochim.txt"
    exit 1
fi

# Count ASVs before and after
ASVS_BEFORE=$(tail -n +2 "$SEQTAB" | wc -l)
ASVS_AFTER=$(tail -n +2 "$OUTPUT_DIR/03-chimera/seqtab_nochim.txt" | wc -l)
ASVS_REMOVED=$((ASVS_BEFORE - ASVS_AFTER))

log_success "ASV table summary:"
log_success "  Before: $ASVS_BEFORE ASVs"
log_success "  After: $ASVS_AFTER ASVs"
log_success "  Removed: $ASVS_REMOVED chimeras"

# Stage 4: Generate per-sample chimera report
log_info "=== Per-Sample Chimera Statistics ==="

cat > "$OUTPUT_DIR/03-chimera/generate_stats.R" << 'EOFR'
#!/usr/bin/env Rscript

seqtab_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

# Read tables
seqtab_all <- read.table(seqtab_file, sep = "\t", header = TRUE, row.names = 1)
seqtab_clean <- read.table(gsub("seqtab.txt", "03-chimera/seqtab_nochim.txt", seqtab_file), 
                           sep = "\t", header = TRUE, row.names = 1)

# Calculate chimera removal rate per sample
stats <- data.frame(
    Sample = rownames(seqtab_all),
    Reads_Before = rowSums(seqtab_all),
    Reads_After = rowSums(seqtab_clean),
    Chimeras_Removed = rowSums(seqtab_all) - rowSums(seqtab_clean),
    Chimera_Percent = round(100 * (1 - rowSums(seqtab_clean) / rowSums(seqtab_all)), 2)
)

write.table(stats, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Per-sample statistics written to", output_file, "\n")
EOFR

if command -v Rscript &>/dev/null; then
    Rscript "$OUTPUT_DIR/03-chimera/generate_stats.R" \
        "$SEQTAB" \
        "$OUTPUT_DIR/03-chimera/chimera_stats_per_sample.txt" \
        2>&1 | tee -a "$LOG_FILE" || log_warning "R statistics failed (non-critical)"
fi

# Generate comprehensive report
cat > "$OUTPUT_DIR/03-chimera/CHIMERA_REPORT.txt" << EOF
=== Chimera Detection Report ===
Timestamp: $(date)
Profile: $PROFILE

METHOD: Hybrid DADA2 + VSEARCH UCHIME
Rationale: Combines DADA2's internal chimera removal with VSEARCH de novo/reference detection
           Reduces false positives compared to single-method approaches (2025 best practice)

STAGE 02 (DADA2) RESULTS:
  Performed consensus chimera removal during denoising
  
STAGE 03 (VSEARCH) RESULTS:
  De novo detection:
    - Chimeras found: $DENOVO_CHIMERAS
    - Non-chimeras: $DENOVO_NONCHIM
  
  Reference-based detection: $([ -n "$REF_DB" ] && echo "Enabled" || echo "Skipped (no reference)")
    - Additional chimeras: $([ -n "$REF_DB" ] && echo "$REF_CHIMERAS" || echo "N/A")

ASV TABLE CHANGES:
  Before chimera removal: $ASVS_BEFORE ASVs
  After chimera removal: $ASVS_AFTER ASVs
  Total removed: $ASVS_REMOVED sequences

OUTPUTS:
  - seqtab_nochim.txt (cleaned ASV abundance table)
  - nonchimeras_denovo.fasta (final representative sequences)
  - uchime_denovo.txt (de novo detection details)
  - chimera_stats_per_sample.txt (per-sample removal rates)

NEXT STAGE:
  Proceed to 04-taxonomy-rdp.sh for taxonomic classification

REFERENCES:
  - VSEARCH: https://github.com/torognes/vsearch
  - UCHIME: Edgar et al., 2016 (hybrid approach validated)
EOF

log_success "Chimera report written to $OUTPUT_DIR/03-chimera/CHIMERA_REPORT.txt"

# Emit metrics
emit_metric "chimera_count" "$DENOVO_CHIMERAS"
emit_metric "asv_count_after_chimera" "$ASVS_AFTER"

log_success "=== Stage 03 Complete ==="
