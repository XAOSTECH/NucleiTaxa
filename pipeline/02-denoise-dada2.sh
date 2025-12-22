#!/bin/bash

###############################################################################
# NucleiTaxa Pipeline Stage 02: Denoise with DADA2
# 
# Amplicon Sequence Variant (ASV) inference using DADA2's error correction
# DADA2 (2023 best practice for 16S/ITS) outperforms OTU clustering
#
# Input:  Filtered paired-end FASTQ files from Stage 01
# Output: ASV table, representative sequences, error model plots
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

# DADA2 parameters (configurable)
TRUNCATE_LENGTH_F=240
TRUNCATE_LENGTH_R=160
TRIM_LEFT_F=0
TRIM_LEFT_R=0
MAX_EE_F=2
MAX_EE_R=2
MIN_OVERLAP=12

# Reference database for taxonomy (optional, for reference-based QC)
REF_FASTA=""

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
        --input-dir)    INPUT_DIR="$2"; shift 2 ;;
        --output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
        --profile)      PROFILE="$2"; shift 2 ;;
        --jobs)         JOBS="$2"; shift 2 ;;
        --config)       CONFIG_FILE="$2"; shift 2 ;;
        --analytics-port) ANALYTICS_PORT="$2"; shift 2 ;;
        *)              log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Usage: $0 --input-dir <path> --output-dir <path>"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/02-denoise"
LOG_FILE="$OUTPUT_DIR/02-denoise/denoise.log"

log_info "=== NucleiTaxa Stage 02: DADA2 Denoise ==="
log_info "Input:  $INPUT_DIR/01-preprocess"
log_info "Output: $OUTPUT_DIR/02-denoise"
log_info "Profile: $PROFILE"

# Verify DADA2 R package is installed
if ! Rscript -e "library(dada2)" 2>/dev/null; then
    log_error "DADA2 R package not found. Install with: R -e 'install.packages(\"dada2\", repos=\"http://bioconductor.org\")'"
    exit 1
fi

log_success "DADA2 R package verified"

# Load profile settings
if [[ -f "$CONFIG_FILE" ]]; then
    log_info "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE" || true
fi

# Generate R script for DADA2 analysis
DADA2_SCRIPT="$OUTPUT_DIR/02-denoise/dada2_pipeline.R"

cat > "$DADA2_SCRIPT" << 'EOFR'
#!/usr/bin/env Rscript

# DADA2 amplicon analysis pipeline
library(dada2)
library(foreach)
library(doParallel)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
n_threads <- as.numeric(args[3])
trunc_f <- as.numeric(args[4])
trunc_r <- as.numeric(args[5])
trim_left_f <- as.numeric(args[6])
trim_left_r <- as.numeric(args[7])
max_ee_f <- as.numeric(args[8])
max_ee_r <- as.numeric(args[9])
min_overlap <- as.numeric(args[10])

cat("=== DADA2 Pipeline ===\n")
cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Threads:", n_threads, "\n")
cat("Truncation lengths (F/R):", trunc_f, "/", trunc_r, "\n")
cat("Quality threshold (maxEE):", max_ee_f, "/", max_ee_r, "\n\n")

# Set up parallel processing
registerDoParallel(cores = n_threads)

# List filtered FASTQ files
fnFs <- sort(list.files(file.path(input_dir, "01-preprocess"), 
                        pattern = "_R1_filtered.fastq.gz", 
                        full.names = TRUE))
fnRs <- sort(list.files(file.path(input_dir, "01-preprocess"), 
                        pattern = "_R2_filtered.fastq.gz", 
                        full.names = TRUE))

if (length(fnFs) == 0) {
    stop("No R1 filtered FASTQ files found!")
}

cat("Found", length(fnFs), "paired-end samples\n\n")

# Quality plots (before filtering)
cat("Generating quality plots...\n")
pdf(file.path(output_dir, "quality_plots_before.pdf"), width = 12, height = 10)
plotQualityProfile(fnFs[1:min(4, length(fnFs))])
dev.off()
cat("Quality plot saved\n")

# Filter and trim
cat("\n=== Filtering and Trimming ===\n")
filtFs <- file.path(output_dir, "filtered", basename(fnFs))
filtRs <- file.path(output_dir, "filtered", basename(fnRs))
dir.create(file.path(output_dir, "filtered"), showWarnings = FALSE)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(trunc_f, trunc_r),
                     trimLeft = c(trim_left_f, trim_left_r),
                     maxEE = c(max_ee_f, max_ee_r),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = n_threads,
                     verbose = TRUE)

cat("\nFilterAndTrim results:\n")
print(out)

# Error model learning
cat("\n=== Learning Error Rates ===\n")
errF <- learnErrors(filtFs, multithread = n_threads, verbose = TRUE)
errR <- learnErrors(filtRs, multithread = n_threads, verbose = TRUE)

pdf(file.path(output_dir, "error_plots.pdf"), width = 12, height = 10)
plotErrors(errF, ask = FALSE, main = "Forward Error Rates")
plotErrors(errR, ask = FALSE, main = "Reverse Error Rates")
dev.off()
cat("Error plots saved\n")

# Dereplicate sequences
cat("\n=== Dereplication ===\n")
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sapply(strsplit(basename(fnFs), "_R1"), "[", 1)
names(derepRs) <- sapply(strsplit(basename(fnRs), "_R2"), "[", 1)

# DADA2 sample inference
cat("\n=== Sample Inference (DADA2) ===\n")
dadaFs <- dada(derepFs, err = errF, multithread = n_threads, verbose = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = n_threads, verbose = TRUE)

cat("DADA2 inference complete\n")
cat("Forward ASVs:", sum(sapply(dadaFs, function(x) nrow(x$clustering))), "\n")
cat("Reverse ASVs:", sum(sapply(dadaRs, function(x) nrow(x$clustering))), "\n")

# Merge paired reads
cat("\n=== Merging Paired Reads ===\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                      minOverlap = min_overlap,
                      verbose = TRUE,
                      propagateCol = c("n0", "n1", "n2", "n3", "n4", "n5"))

# Create ASV table
cat("\n=== Creating ASV Table ===\n")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
cat("ASV table dimensions:", nrow(seqtab), "samples x", ncol(seqtab), "ASVs\n")

# Remove chimeras
cat("\n=== Chimera Removal (DADA2) ===\n")
seqtab.nochim <- removeBimeras(seqtab, method = "consensus", multithread = n_threads, verbose = TRUE)
cat("Chimeras identified:", ncol(seqtab) - ncol(seqtab.nochim), "\n")
cat("Final ASV count:", ncol(seqtab.nochim), "\n")

# Quality summary
cat("\n=== Read Retention Summary ===\n")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN),
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- names(derepFs)
head(track)
write.table(track, file.path(output_dir, "read_retention_track.txt"), sep = "\t", quote = FALSE)

# Save outputs
cat("\n=== Saving Outputs ===\n")
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab.rds"))
write.table(seqtab.nochim, file.path(output_dir, "seqtab.txt"), sep = "\t", quote = FALSE)

# Export representative sequences (FASTA)
uniquesToFasta(seqtab.nochim, file.path(output_dir, "asv_sequences.fasta"))
cat("ASV sequences saved to asv_sequences.fasta\n")

# Export taxonomy-ready format
cat("\nPipeline complete!\n")
cat("Outputs in", output_dir, ":\n")
cat("  - seqtab.txt (ASV abundance table)\n")
cat("  - asv_sequences.fasta (representative sequences)\n")
cat("  - read_retention_track.txt (quality summary)\n")
cat("  - error_plots.pdf (error model validation)\n")
EOFR

log_success "DADA2 R script generated"

# Execute DADA2 pipeline
log_info "Running DADA2 pipeline..."

Rscript "$DADA2_SCRIPT" \
    "$OUTPUT_DIR" \
    "$OUTPUT_DIR/02-denoise" \
    "$JOBS" \
    "$TRUNCATE_LENGTH_F" \
    "$TRUNCATE_LENGTH_R" \
    "$TRIM_LEFT_F" \
    "$TRIM_LEFT_R" \
    "$MAX_EE_F" \
    "$MAX_EE_R" \
    "$MIN_OVERLAP" \
    2>&1 | tee -a "$LOG_FILE"

if [[ ! -f "$OUTPUT_DIR/02-denoise/seqtab.txt" ]]; then
    log_error "DADA2 failed: seqtab.txt not generated"
    exit 1
fi

# Count ASVs
ASV_COUNT=$(tail -n +2 "$OUTPUT_DIR/02-denoise/seqtab.txt" | wc -l)
log_success "DADA2 complete: $ASV_COUNT ASVs inferred"

emit_metric "asv_count" "$ASV_COUNT"

# Generate report
cat > "$OUTPUT_DIR/02-denoise/DADA2_REPORT.txt" << EOF
=== DADA2 Denoising Report ===
Timestamp: $(date)
Profile: $PROFILE

PARAMETERS:
  Truncation (F/R): ${TRUNCATE_LENGTH_F}/${TRUNCATE_LENGTH_R} bp
  Trimming (F/R): ${TRIM_LEFT_F}/${TRIM_LEFT_R} bp
  Max expected errors (F/R): ${MAX_EE_F}/${MAX_EE_R}
  Min overlap: $MIN_OVERLAP bp

RESULTS:
  Total ASVs: $ASV_COUNT
  Output files:
    - seqtab.txt (abundance table)
    - asv_sequences.fasta (representative sequences)
    - read_retention_track.txt (per-sample statistics)
    - error_plots.pdf (model validation)
    - quality_plots_before.pdf (input quality)

NEXT STAGE:
  Proceed to 03-chimera-vsearch.sh for additional chimera detection
EOF

log_success "DADA2 report written"
log_success "=== Stage 02 Complete ==="
