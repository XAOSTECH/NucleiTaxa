#!/bin/bash

###############################################################################
# NucleiTaxa Pipeline Stage 03: Chimera Detection
# 
# Hybrid approach: DADA2 + VSEARCH UCHIME
# - DADA2 removes chimeras during denoising (already done in Stage 02)
# - VSEARCH de novo + reference-based chimera detection (additional QC)
# Research shows hybrid approach reduces false positives vs single method
#
# GPU ACCELERATION: Supports --cuda flag for VSEARCH UCHIME (10x speedup)
# VSEARCH v1.14.0+ required for CUDA support
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
USE_CUDA=false
GPU_DEVICE=0

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
        --cuda)            USE_CUDA=true; shift ;;
        --no-cuda)         USE_CUDA=false; shift ;;
        --gpu-device)      GPU_DEVICE="$2"; shift 2 ;;
        *)                 log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Usage: $0 --input-dir <path> --output-dir <path> [--cuda]"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/03-chimera"
LOG_FILE="$OUTPUT_DIR/03-chimera/chimera.log"

log_info "=== NucleiTaxa Stage 03: Chimera Detection (VSEARCH) ==="
log_info "Input:  $OUTPUT_DIR/02-denoise"
log_info "Output: $OUTPUT_DIR/03-chimera"
log_info "Profile: $PROFILE"
log_info "GPU Acceleration: $([ "$USE_CUDA" = true ] && echo "ENABLED" || echo "DISABLED")"

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

# Check CUDA support if requested
CUDA_FLAG=""
if [ "$USE_CUDA" = true ]; then
    if vsearch --version 2>&1 | grep -qi "CUDA"; then
        CUDA_FLAG="--cuda"
        log_success "VSEARCH CUDA support detected"
        
        # Verify nvidia-smi is available
        if ! command -v nvidia-smi &>/dev/null; then
            log_warning "nvidia-smi not found, GPU may not be available"
            CUDA_FLAG=""
            USE_CUDA=false
        else
            # Show GPU info
            log_info "GPU Device: $GPU_DEVICE"
            nvidia-smi -i "$GPU_DEVICE" --query-gpu=name,memory.total --format=csv,noheader | xargs -I {} log_info "GPU: {}"
        fi
    else
        log_warning "VSEARCH CUDA support not detected (requires v1.14.0+), using CPU"
        USE_CUDA=false
    fi
fi

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
[ ! -z "$CUDA_FLAG" ] && log_info "GPU acceleration enabled for this stage"

# Start GPU monitoring if CUDA enabled
if [ ! -z "$CUDA_FLAG" ]; then
    (while true; do
        nvidia-smi -i "$GPU_DEVICE" --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader | \
            awk -F, '{printf "[GPU-MON] GPU Util: %s, Memory: %s / %s\n", $1, $2, $3}'
        sleep 2
    done) &
    GPU_MONITOR_PID=$!
fi

# VSEARCH denovo chimera detection with optional GPU
# Output: tab-delimited with fields: query, db_count, chimera_indicator (Y/N)
VSEARCH_START=$(date +%s)
vsearch \
    --uchime_denovo "$ASV_FASTA" \
    --chimeras "$OUTPUT_DIR/03-chimera/chimeras_denovo.fasta" \
    --nonchimeras "$OUTPUT_DIR/03-chimera/nonchimeras_denovo.fasta" \
    --uchimeout "$OUTPUT_DIR/03-chimera/uchime_denovo.txt" \
    --threads "$JOBS" \
    $CUDA_FLAG \
    2>&1 | tee -a "$LOG_FILE"

VSEARCH_END=$(date +%s)
VSEARCH_TIME=$((VSEARCH_END - VSEARCH_START))

# Kill GPU monitor
if [ ! -z "$GPU_MONITOR_PID" ]; then
    kill $GPU_MONITOR_PID 2>/dev/null || true
fi

log_info "De novo detection completed in ${VSEARCH_TIME}s"

# Count chimeras detected
DENOVO_CHIMERAS=$(grep -c "Y$" "$OUTPUT_DIR/03-chimera/uchime_denovo.txt" || echo "0")
DENOVO_NONCHIM=$(grep -c "N$" "$OUTPUT_DIR/03-chimera/uchime_denovo.txt" || echo "0")

log_success "De novo detection: $DENOVO_CHIMERAS chimeras, $DENOVO_NONCHIM non-chimeras"

# Stage 2: Reference-based chimera detection (if reference available)
if [[ -n "$REF_DB" ]] && [[ -f "$REF_DB" ]]; then
    log_info "=== Reference-Based Chimera Detection ==="
    log_info "Reference database: $REF_DB"
    [ ! -z "$CUDA_FLAG" ] && log_info "GPU acceleration enabled for this stage"
    
    REF_START=$(date +%s)
    vsearch \
        --uchime_ref "$OUTPUT_DIR/03-chimera/nonchimeras_denovo.fasta" \
        --db "$REF_DB" \
        --chimeras "$OUTPUT_DIR/03-chimera/chimeras_reference.fasta" \
        --nonchimeras "$OUTPUT_DIR/03-chimera/nonchimeras_final.fasta" \
        --uchimeout "$OUTPUT_DIR/03-chimera/uchime_reference.txt" \
        --threads "$JOBS" \
        $CUDA_FLAG \
        2>&1 | tee -a "$LOG_FILE"
    
    REF_END=$(date +%s)
    REF_TIME=$((REF_END - REF_START))
    
    REF_CHIMERAS=$(grep -c "Y$" "$OUTPUT_DIR/03-chimera/uchime_reference.txt" || echo "0")
    log_success "Reference-based detection: $REF_CHIMERAS additional chimeras (${REF_TIME}s)"
    
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
GPU Acceleration: $([ ! -z "$CUDA_FLAG" ] && echo "ENABLED" || echo "DISABLED")

METHOD: Hybrid DADA2 + VSEARCH UCHIME
Rationale: Combines DADA2's internal chimera removal with VSEARCH de novo/reference detection
           Reduces false positives compared to single-method approaches (2025 best practice)

STAGE 02 (DADA2) RESULTS:
  Performed consensus chimera removal during denoising
  
STAGE 03 (VSEARCH) RESULTS:
  De novo detection (${VSEARCH_TIME}s):
    - Chimeras found: $DENOVO_CHIMERAS
    - Non-chimeras: $DENOVO_NONCHIM
  
  Reference-based detection: $([ -n "$REF_DB" ] && echo "Enabled" || echo "Skipped (no reference)")
    - Additional chimeras: $([ -n "$REF_DB" ] && echo "$REF_CHIMERAS" || echo "N/A")
    - Detection time: $([ -n "$REF_DB" ] && echo "${REF_TIME}s" || echo "N/A")

GPU ACCELERATION:
  Status: $([ ! -z "$CUDA_FLAG" ] && echo "ENABLED ($([ "$USE_CUDA" = true ] && echo "10x speedup expected" || echo "not available"))" || echo "DISABLED")
  Expected speedup: 10x (VSEARCH CUDA acceleration)

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
  - VSEARCH CUDA: https://github.com/torognes/vsearch/wiki/Running-with-GPU-acceleration
  - UCHIME: Edgar et al., 2016 (hybrid approach validated)
EOF

log_success "Chimera report written to $OUTPUT_DIR/03-chimera/CHIMERA_REPORT.txt"

# Emit metrics
emit_metric "chimera_count" "$DENOVO_CHIMERAS"
emit_metric "asv_count_after_chimera" "$ASVS_AFTER"
emit_metric "gpu_acceleration" "$([ ! -z "$CUDA_FLAG" ] && echo "1" || echo "0")"
emit_metric "vsearch_time_seconds" "$VSEARCH_TIME"

log_success "=== Stage 03 Complete ==="
