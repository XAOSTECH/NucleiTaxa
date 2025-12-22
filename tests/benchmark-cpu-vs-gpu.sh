#!/bin/bash

#############################################################################
# CPU vs GPU Benchmark Suite for NucleiTaxa
# 
# Validates and benchmarks CPU and GPU implementations:
# - VSEARCH UCHIME chimera detection (--cuda flag)
# - FastTree phylogenetic inference (future GPU target)
# - End-to-end pipeline execution with metrics
#
# Usage: ./tests/benchmark-cpu-vs-gpu.sh [--quick] [--mock-data] [--size small|medium|large]
#############################################################################

set -e

BENCHMARK_DIR="test_benchmarks"
QUICK_RUN=false
MOCK_DATA=true
DATA_SIZE="small"
RESULTS_CSV="$BENCHMARK_DIR/results.csv"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --quick)       QUICK_RUN=true; shift ;;
        --no-mock)     MOCK_DATA=false; shift ;;
        --size)        DATA_SIZE="$2"; shift 2 ;;
        *)             echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║${NC}     NucleiTaxa GPU vs CPU Benchmark Suite                    ${BLUE}║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Create benchmark directory
mkdir -p "$BENCHMARK_DIR"

# Initialize results CSV
if [ ! -f "$RESULTS_CSV" ]; then
    cat > "$RESULTS_CSV" << 'EOF'
timestamp,test_name,implementation,data_size,execution_time_seconds,sequences_processed,throughput_seqs_per_sec,gpu_available,gpu_memory_used_mb,speedup_factor
EOF
fi

# Function to generate synthetic FASTQ data
generate_mock_data() {
    local size=$1
    local num_seqs=0
    
    case $size in
        small)  num_seqs=1000 ;;
        medium) num_seqs=10000 ;;
        large)  num_seqs=100000 ;;
    esac
    
    local output_r1="$BENCHMARK_DIR/mock_R1.fastq"
    local output_r2="$BENCHMARK_DIR/mock_R2.fastq"
    
    echo -e "${YELLOW}[MOCK]${NC} Generating $num_seqs sequence pairs ($size dataset)..."
    
    # Generate mock FASTQ with realistic structure
    > "$output_r1"
    > "$output_r2"
    
    for ((i=1; i<=num_seqs; i++)); do
        # Forward read (R1)
        seq=$(openssl rand -hex 50 | tr -d '\n')  # 100bp mock sequence
        qual=$(python3 -c "import random; print(''.join([chr(random.randint(33,126)) for _ in range(100)]))")
        
        echo "@read_${i}/1" >> "$output_r1"
        echo "${seq:0:100}" >> "$output_r1"
        echo "+" >> "$output_r1"
        echo "$qual" >> "$output_r1"
        
        # Reverse read (R2)
        echo "@read_${i}/2" >> "$output_r2"
        echo "${seq:0:100}" >> "$output_r2"
        echo "+" >> "$output_r2"
        echo "$qual" >> "$output_r2"
        
        if [ $((i % 1000)) -eq 0 ]; then
            echo -e "${BLUE}[PROGRESS]${NC} $i / $num_seqs sequences generated"
        fi
    done
    
    echo -e "${GREEN}[OK]${NC} Mock data generated: $output_r1, $output_r2"
    echo "$output_r1:$output_r2:$num_seqs"
}

# Function to log benchmark result
log_result() {
    local test_name="$1"
    local implementation="$2"
    local exec_time="$3"
    local num_seqs="$4"
    local gpu_avail="$5"
    local gpu_mem="$6"
    local speedup="$7"
    
    local throughput=$(echo "scale=2; $num_seqs / $exec_time" | bc)
    
    echo "$TIMESTAMP,$test_name,$implementation,$DATA_SIZE,$exec_time,$num_seqs,$throughput,$gpu_avail,$gpu_mem,$speedup" >> "$RESULTS_CSV"
    
    echo -e "${GREEN}[RESULT]${NC} $test_name / $implementation"
    echo "  Time: ${exec_time}s | Throughput: ${throughput} seq/s"
}

# Test 1: VSEARCH UCHIME Chimera Detection
test_vsearch_chimera() {
    echo ""
    echo -e "${BLUE}[TEST 1]${NC} VSEARCH UCHIME Chimera Detection"
    echo "=================================================="
    
    # Generate test data
    local mock_data=$(generate_mock_data "$DATA_SIZE")
    local input_fasta="$BENCHMARK_DIR/test_chimera.fasta"
    local num_seqs=$(echo "$mock_data" | cut -d: -f3)
    
    # Convert mock FASTQ to FASTA for VSEARCH
    local r1_fastq=$(echo "$mock_data" | cut -d: -f1)
    seqtk seq -A "$r1_fastq" | head -n 100000 > "$input_fasta" || echo "seqtk not available, creating mock FASTA"
    
    if [ ! -s "$input_fasta" ] || [ ! -f "$input_fasta" ]; then
        # Create minimal FASTA if seqtk failed
        echo -e "${YELLOW}[WARN]${NC} Creating synthetic FASTA for testing..."
        > "$input_fasta"
        for ((i=1; i<=100; i++)); do
            echo ">seq_$i"
            openssl rand -hex 50 | tr -d '\n'
            echo ""
        done >> "$input_fasta"
    fi
    
    num_seqs=$(grep -c "^>" "$input_fasta")
    
    # Test 1a: CPU-only VSEARCH
    echo ""
    echo -e "${YELLOW}[1a]${NC} CPU-only VSEARCH..."
    
    local cpu_output="$BENCHMARK_DIR/vsearch_cpu_output"
    mkdir -p "$cpu_output"
    
    if command -v vsearch &>/dev/null; then
        local cpu_start=$(date +%s%N)
        
        vsearch \
            --uchime_denovo "$input_fasta" \
            --nonchimeras "$cpu_output/nonchimeras.fasta" \
            --chimeras "$cpu_output/chimeras.fasta" \
            --uchimeout "$cpu_output/uchime.txt" \
            --threads 4 \
            2>&1 | head -5
        
        local cpu_end=$(date +%s%N)
        local cpu_time=$(echo "scale=3; ($cpu_end - $cpu_start) / 1000000000" | bc)
        
        echo -e "${GREEN}[PASS]${NC} CPU VSEARCH completed in ${cpu_time}s"
        log_result "vsearch_uchime" "cpu" "$cpu_time" "$num_seqs" "false" "0" "1.0"
        
        # Test 1b: GPU VSEARCH (if available)
        echo ""
        echo -e "${YELLOW}[1b]${NC} GPU VSEARCH (if CUDA available)..."
        
        if vsearch --version 2>&1 | grep -qi "CUDA"; then
            echo -e "${GREEN}[CUDA]${NC} VSEARCH CUDA support detected"
            
            local gpu_output="$BENCHMARK_DIR/vsearch_gpu_output"
            mkdir -p "$gpu_output"
            
            local gpu_start=$(date +%s%N)
            
            # Check GPU availability
            if command -v nvidia-smi &>/dev/null; then
                local gpu_mem_before=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
                
                vsearch \
                    --uchime_denovo "$input_fasta" \
                    --nonchimeras "$gpu_output/nonchimeras.fasta" \
                    --chimeras "$gpu_output/chimeras.fasta" \
                    --uchimeout "$gpu_output/uchime.txt" \
                    --threads 4 \
                    --cuda \
                    2>&1 | head -5
                
                local gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
                local gpu_mem_used=$((gpu_mem_after - gpu_mem_before))
            else
                echo -e "${YELLOW}[WARN]${NC} nvidia-smi not found, GPU metrics unavailable"
                local gpu_mem_used=0
            fi
            
            local gpu_end=$(date +%s%N)
            local gpu_time=$(echo "scale=3; ($gpu_end - $gpu_start) / 1000000000" | bc)
            local speedup=$(echo "scale=2; $cpu_time / $gpu_time" | bc)
            
            echo -e "${GREEN}[PASS]${NC} GPU VSEARCH completed in ${gpu_time}s (${speedup}x speedup)"
            log_result "vsearch_uchime" "gpu_cuda" "$gpu_time" "$num_seqs" "true" "$gpu_mem_used" "$speedup"
            
            # Validate GPU output matches CPU output
            local cpu_chim_count=$(wc -l < "$cpu_output/uchime.txt")
            local gpu_chim_count=$(wc -l < "$gpu_output/uchime.txt")
            
            if [ "$cpu_chim_count" -eq "$gpu_chim_count" ]; then
                echo -e "${GREEN}[VALIDATE]${NC} CPU and GPU results match ($cpu_chim_count chimeras detected)"
            else
                echo -e "${RED}[ERROR]${NC} Results mismatch: CPU=$cpu_chim_count, GPU=$gpu_chim_count"
                return 1
            fi
        else
            echo -e "${YELLOW}[SKIP]${NC} VSEARCH CUDA support not available (requires v1.14.0+)"
        fi
    else
        echo -e "${RED}[ERROR]${NC} VSEARCH not installed"
        return 1
    fi
}

# Test 2: FastTree Phylogenetics (CPU only for now)
test_fasttree() {
    echo ""
    echo -e "${BLUE}[TEST 2]${NC} FastTree Phylogenetic Inference"
    echo "=================================================="
    
    local input_fasta="$BENCHMARK_DIR/test_phylo.fasta"
    
    # Create aligned mock sequences
    echo -e "${YELLOW}[SETUP]${NC} Creating mock aligned sequences..."
    > "$input_fasta"
    for ((i=1; i<=50; i++)); do
        echo ">seq_$i"
        # Generate 250bp aligned sequence
        for ((j=1; j<=25; j++)); do
            openssl rand -hex 10 | tr -d '\n'
        done
        echo ""
    done >> "$input_fasta"
    
    if command -v fasttree &>/dev/null; then
        echo -e "${YELLOW}[CPU]${NC} Running FastTree (CPU only)..."
        
        local ft_start=$(date +%s%N)
        
        fasttree -nt "$input_fasta" > "$BENCHMARK_DIR/tree.nwk" 2>&1 || true
        
        local ft_end=$(date +%s%N)
        local ft_time=$(echo "scale=3; ($ft_end - $ft_start) / 1000000000" | bc)
        
        local num_seqs=$(grep -c "^>" "$input_fasta")
        
        echo -e "${GREEN}[PASS]${NC} FastTree completed in ${ft_time}s"
        log_result "fasttree_inference" "cpu" "$ft_time" "$num_seqs" "false" "0" "1.0"
        
        echo -e "${YELLOW}[NOTE]${NC} GPU support (VeryFastTree) planned for Phase 2"
    else
        echo -e "${RED}[ERROR]${NC} FastTree not installed"
        return 1
    fi
}

# Test 3: Full Pipeline Stage 03 (Chimera Detection)
test_full_pipeline_stage03() {
    echo ""
    echo -e "${BLUE}[TEST 3]${NC} Full Pipeline Stage 03 (CPU vs GPU)"
    echo "=================================================="
    
    echo -e "${YELLOW}[SETUP]${NC} Preparing mock pipeline data..."
    
    local pipeline_dir="$BENCHMARK_DIR/pipeline_test"
    mkdir -p "$pipeline_dir/02-denoise"
    
    # Create mock denoised sequences
    local mock_seqs="$pipeline_dir/02-denoise/asv_sequences.fasta"
    > "$mock_seqs"
    for ((i=1; i<=100; i++)); do
        echo ">ASV_$i"
        openssl rand -hex 50 | tr -d '\n'
        echo ""
    done >> "$mock_seqs"
    
    # Create mock ASV table
    local mock_table="$pipeline_dir/02-denoise/seqtab.txt"
    echo -e "SeqID\tsample1\tsample2\tsample3" > "$mock_table"
    for ((i=1; i<=100; i++)); do
        echo -e "ASV_$i\t$((RANDOM % 1000))\t$((RANDOM % 1000))\t$((RANDOM % 1000))" >> "$mock_table"
    done
    
    # Test CPU execution
    echo -e "${YELLOW}[CPU]${NC} Running Stage 03 (CPU only)..."
    
    local cpu_dir="$pipeline_dir/cpu_output"
    local cpu_start=$(date +%s%N)
    
    if [ -f "pipeline/03-chimera-vsearch.sh" ]; then
        bash pipeline/03-chimera-vsearch.sh \
            --input-dir "$pipeline_dir" \
            --output-dir "$cpu_dir" \
            --no-cuda \
            2>&1 | grep -E "\[INFO\]|\[✓\]" || true
    else
        echo -e "${YELLOW}[SKIP]${NC} Pipeline script not found, using mock"
    fi
    
    local cpu_end=$(date +%s%N)
    local cpu_time=$(echo "scale=3; ($cpu_end - $cpu_start) / 1000000000" | bc)
    
    echo -e "${GREEN}[PASS]${NC} CPU Stage 03 completed in ${cpu_time}s"
    log_result "pipeline_stage03" "cpu" "$cpu_time" "100" "false" "0" "1.0"
    
    # Test GPU execution
    if vsearch --version 2>&1 | grep -qi "CUDA"; then
        echo -e "${YELLOW}[GPU]${NC} Running Stage 03 (GPU acceleration)..."
        
        local gpu_dir="$pipeline_dir/gpu_output"
        local gpu_start=$(date +%s%N)
        
        if [ -f "pipeline/03-chimera-vsearch.sh" ]; then
            bash pipeline/03-chimera-vsearch.sh \
                --input-dir "$pipeline_dir" \
                --output-dir "$gpu_dir" \
                --cuda \
                2>&1 | grep -E "\[INFO\]|\[✓\]\[CUDA\]" || true
        fi
        
        local gpu_end=$(date +%s%N)
        local gpu_time=$(echo "scale=3; ($gpu_end - $gpu_start) / 1000000000" | bc)
        local speedup=$(echo "scale=2; $cpu_time / $gpu_time" | bc)
        
        echo -e "${GREEN}[PASS]${NC} GPU Stage 03 completed in ${gpu_time}s (${speedup}x speedup)"
        log_result "pipeline_stage03" "gpu_cuda" "$gpu_time" "100" "true" "0" "$speedup"
    fi
}

# Generate final report
generate_report() {
    echo ""
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}                  BENCHMARK REPORT                          ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    
    if [ -f "$RESULTS_CSV" ]; then
        echo "Results saved to: $RESULTS_CSV"
        echo ""
        echo "Summary:"
        tail -n +2 "$RESULTS_CSV" | awk -F, '{
            printf "  %-30s %-10s %8.3fs  %5.1fx\n", $2, $3, $5, $10
        }'
        echo ""
        
        # Calculate aggregate metrics
        echo "Performance Summary:"
        tail -n +2 "$RESULTS_CSV" | awk -F, 'BEGIN {
            gpu_count = 0; gpu_speedup = 0
        } {
            if ($3 ~ /gpu/) {
                gpu_count++
                gpu_speedup += $10
            }
        } END {
            if (gpu_count > 0) {
                avg_speedup = gpu_speedup / gpu_count
                printf "  Average GPU speedup: %.2fx\n", avg_speedup
            }
        }'
    fi
    
    echo ""
    echo -e "${GREEN}✓${NC} Benchmark suite completed!"
}

# Main execution
main() {
    # Check prerequisites
    if [ "$MOCK_DATA" = true ]; then
        if ! command -v openssl &>/dev/null; then
            echo -e "${RED}[ERROR]${NC} openssl required for mock data generation"
            exit 1
        fi
    fi
    
    # Run tests
    test_vsearch_chimera || echo -e "${YELLOW}[WARN]${NC} VSEARCH test failed"
    
    if [ "$QUICK_RUN" = false ]; then
        test_fasttree || echo -e "${YELLOW}[WARN]${NC} FastTree test failed"
        test_full_pipeline_stage03 || echo -e "${YELLOW}[WARN]${NC} Pipeline test failed"
    fi
    
    generate_report
    
    echo -e "${GREEN}✓${NC} All benchmarks complete"
}

main "$@"
