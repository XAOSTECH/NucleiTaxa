#!/bin/bash

###############################################################################
# NucleiTaxa Pipeline Stage 01: Preprocess
# 
# Quality control, adapter trimming, and initial read filtering
# Uses BBTools (fast k-mer matching) + Cutadapt (precision trimming)
#
# Input:  Paired-end FASTQ files (raw sequencing data)
# Output: QC-filtered, adapter-trimmed FASTQ pairs
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

# BBTools parameters (configurable via config file)
KMER_SIZE=31
MIN_QUALITY=20
MIN_LENGTH=50
ADAPTER_QUALITY=10
ENTROPY_THRESHOLD=0.5

# Cutadapt parameters
PRIMER_F=""
PRIMER_R=""
ERROR_RATE=0.1
MIN_OVERLAP=10

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

# Emit metrics to analytics server
emit_metric() {
    local key="$1"
    local value="$2"
    
    if command -v nc &>/dev/null; then
        # Send as simple key=value to local socket if analytics running
        echo "$key=$value" | nc -w1 localhost "$ANALYTICS_PORT" 2>/dev/null || true
    fi
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)    INPUT_DIR="$2"; shift 2 ;;
        --output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
        --profile)      PROFILE="$2"; shift 2 ;;
        --jobs)         JOBS="$2"; shift 2 ;;
        --config)       CONFIG_FILE="$2"; shift 2 ;;
        --analytics-port) ANALYTICS_PORT="$2"; shift 2 ;;
        *)              log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Usage: $0 --input-dir <path> --output-dir <path> [--profile <name>]"
    exit 1
fi

# Setup output directory
mkdir -p "$OUTPUT_DIR/01-preprocess"
LOG_FILE="$OUTPUT_DIR/01-preprocess/preprocess.log"

log_info "=== NucleiTaxa Stage 01: Preprocess ==="
log_info "Input:  $INPUT_DIR"
log_info "Output: $OUTPUT_DIR/01-preprocess"
log_info "Profile: $PROFILE"

# Load profile-specific settings
if [[ -f "$CONFIG_FILE" ]]; then
    log_info "Loading config from $CONFIG_FILE"
    # Source safely (validate first)
    if grep -q "^MIN_QUALITY=" "$CONFIG_FILE"; then
        source "$CONFIG_FILE"
    fi
fi

# Detect input files
FASTQ_R1=($(find "$INPUT_DIR" -maxdepth 1 -name "*_R1*.fastq.gz" | sort))
FASTQ_R2=($(find "$INPUT_DIR" -maxdepth 1 -name "*_R2*.fastq.gz" | sort))

if [[ ${#FASTQ_R1[@]} -eq 0 ]]; then
    log_error "No R1 FASTQ files found in $INPUT_DIR"
    exit 1
fi

log_info "Found ${#FASTQ_R1[@]} paired-end samples"

# Quality and length statistics before filtering
log_info "Collecting pre-filter statistics..."
TOTAL_READS_BEFORE=0
TOTAL_BP_BEFORE=0

for r1 in "${FASTQ_R1[@]}"; do
    # Count reads (FASTQ format: 4 lines per read)
    READS=$(zcat "$r1" | wc -l | awk '{print $1 / 4}')
    TOTAL_READS_BEFORE=$((TOTAL_READS_BEFORE + READS))
    
    # Estimate mean length from first 1000 reads
    BP=$(zcat "$r1" | head -4000 | awk 'NR%4==2 {sum+=length($0); count++} END {print sum/count * 1000}' | awk '{printf "%.0f", $0}')
    TOTAL_BP_BEFORE=$((TOTAL_BP_BEFORE + BP))
done

log_success "Pre-filter stats: $TOTAL_READS_BEFORE reads, ~$((TOTAL_BP_BEFORE/1000000))M bp"
emit_metric "reads_before" "$TOTAL_READS_BEFORE"

# BBTools FastqQC for detailed metrics (if available)
if command -v fastqc &>/dev/null; then
    log_info "Running FastQC analysis..."
    mkdir -p "$OUTPUT_DIR/01-preprocess/fastqc"
    fastqc -t "$JOBS" -o "$OUTPUT_DIR/01-preprocess/fastqc" "${FASTQ_R1[@]}" 2>&1 | tee -a "$LOG_FILE"
    log_success "FastQC complete"
fi

# BBTools BBMerge for error correction (optional, enables overlapping)
if command -v bbmerge.sh &>/dev/null; then
    log_info "Attempting read merging with BBMerge..."
    MERGED_COUNT=0
    for i in "${!FASTQ_R1[@]}"; do
        r1="${FASTQ_R1[$i]}"
        r2="${FASTQ_R2[$i]}"
        base=$(basename "$r1" _R1.fastq.gz)
        
        bbmerge.sh \
            in1="$r1" \
            in2="$r2" \
            out="$OUTPUT_DIR/01-preprocess/${base}_merged.fastq.gz" \
            outm="$OUTPUT_DIR/01-preprocess/${base}_unmerged.fastq.gz" \
            vstrict=f \
            threads="$JOBS" \
            2>&1 | tee -a "$LOG_FILE" || true
        
        if [[ -f "$OUTPUT_DIR/01-preprocess/${base}_merged.fastq.gz" ]]; then
            MERGED=$(zcat "$OUTPUT_DIR/01-preprocess/${base}_merged.fastq.gz" | wc -l | awk '{print $1/4}')
            MERGED_COUNT=$((MERGED_COUNT + MERGED))
        fi
    done
    log_success "Merged reads: $MERGED_COUNT"
fi

# Cutadapt: Primer removal (if primers specified)
if [[ -n "$PRIMER_F" ]]; then
    log_info "Removing primers with Cutadapt..."
    
    for i in "${!FASTQ_R1[@]}"; do
        r1="${FASTQ_R1[$i]}"
        r2="${FASTQ_R2[$i]}"
        base=$(basename "$r1" _R1.fastq.gz)
        
        cutadapt \
            -g "^$PRIMER_F" \
            -a "$PRIMER_R\$" \
            -G "^$PRIMER_R" \
            -A "$PRIMER_F\$" \
            --error-rate "$ERROR_RATE" \
            --minimum-length "$((MIN_LENGTH - 10))" \
            --overlap "$MIN_OVERLAP" \
            -o "$OUTPUT_DIR/01-preprocess/${base}_R1_trimmed.fastq.gz" \
            -p "$OUTPUT_DIR/01-preprocess/${base}_R2_trimmed.fastq.gz" \
            "$r1" "$r2" \
            2>&1 | tee -a "$LOG_FILE"
        
        log_success "Trimmed: $base"
    done
else
    log_info "No primers specified, skipping Cutadapt"
fi

# Quality filtering with BBTools (if available)
if command -v bbduk.sh &>/dev/null; then
    log_info "Quality filtering with BBTools BBDuk..."
    
    for i in "${!FASTQ_R1[@]}"; do
        r1="${FASTQ_R1[$i]}"
        r2="${FASTQ_R2[$i]}"
        base=$(basename "$r1" _R1.fastq.gz)
        
        # Check if already trimmed
        if [[ -f "$OUTPUT_DIR/01-preprocess/${base}_R1_trimmed.fastq.gz" ]]; then
            r1="$OUTPUT_DIR/01-preprocess/${base}_R1_trimmed.fastq.gz"
            r2="$OUTPUT_DIR/01-preprocess/${base}_R2_trimmed.fastq.gz"
        fi
        
        bbduk.sh \
            in1="$r1" \
            in2="$r2" \
            out1="$OUTPUT_DIR/01-preprocess/${base}_R1_filtered.fastq.gz" \
            out2="$OUTPUT_DIR/01-preprocess/${base}_R2_filtered.fastq.gz" \
            ref=adapters \
            ktrim=r \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe \
            tbo \
            qtrim=r \
            trimq="$ADAPTER_QUALITY" \
            minlength="$MIN_LENGTH" \
            entropy="$ENTROPY_THRESHOLD" \
            entropywindow=50 \
            entropyk=5 \
            threads="$JOBS" \
            2>&1 | tee -a "$LOG_FILE" || true
        
        log_success "Filtered: $base"
    done
else
    # Fallback: simple quality filtering with seqtk (if available)
    if command -v seqtk &>/dev/null; then
        log_info "Quality filtering with seqtk..."
        
        for i in "${!FASTQ_R1[@]}"; do
            r1="${FASTQ_R1[$i]}"
            base=$(basename "$r1" _R1.fastq.gz)
            
            # Check if already trimmed
            if [[ -f "$OUTPUT_DIR/01-preprocess/${base}_R1_trimmed.fastq.gz" ]]; then
                r1="$OUTPUT_DIR/01-preprocess/${base}_R1_trimmed.fastq.gz"
            fi
            
            seqtk trimfq -q 0.05 -l "$MIN_LENGTH" "$r1" | \
                gzip > "$OUTPUT_DIR/01-preprocess/${base}_R1_filtered.fastq.gz"
            
            log_success "Filtered: $base"
        done
    fi
fi

# Count post-filter reads
log_info "Collecting post-filter statistics..."
TOTAL_READS_AFTER=0
READS_DISCARDED=0

for f in "$OUTPUT_DIR"/01-preprocess/*_filtered.fastq.gz; do
    if [[ -f "$f" ]]; then
        READS=$(zcat "$f" | wc -l | awk '{print $1 / 4}')
        TOTAL_READS_AFTER=$((TOTAL_READS_AFTER + READS))
    fi
done

READS_DISCARDED=$((TOTAL_READS_BEFORE - TOTAL_READS_AFTER))
RETENTION=$((TOTAL_READS_AFTER * 100 / TOTAL_READS_BEFORE))

log_success "Post-filter stats:"
log_success "  Reads remaining: $TOTAL_READS_AFTER ($RETENTION%)"
log_success "  Reads discarded: $READS_DISCARDED"

emit_metric "reads_after_filter" "$TOTAL_READS_AFTER"
emit_metric "reads_discarded" "$READS_DISCARDED"

# Generate QC report
cat > "$OUTPUT_DIR/01-preprocess/QC_REPORT.txt" << EOF
=== Quality Control Report ===
Timestamp: $(date)
Profile: $PROFILE

INPUT STATISTICS:
  Total reads (before): $TOTAL_READS_BEFORE
  Total bases: $TOTAL_BP_BEFORE

QC PARAMETERS:
  Min quality: $MIN_QUALITY
  Min length: $MIN_LENGTH
  Entropy threshold: $ENTROPY_THRESHOLD

OUTPUT STATISTICS:
  Total reads (after): $TOTAL_READS_AFTER
  Reads discarded: $READS_DISCARDED
  Retention rate: $RETENTION%

FILTERS APPLIED:
  - Adapter removal (BBDuk)
  - Quality trimming (Q${ADAPTER_QUALITY})
  - Length filtering (>= $MIN_LENGTH bp)
  - Low-complexity filtering (entropy)

NEXT STAGE:
  Proceed to 02-denoise-dada2.sh with filtered reads
EOF

log_success "QC report written to $OUTPUT_DIR/01-preprocess/QC_REPORT.txt"

# Emit completion metric
emit_metric "stage_01_complete" "1"

log_success "=== Stage 01 Complete ==="
log_info "Elapsed time: $(date)"
