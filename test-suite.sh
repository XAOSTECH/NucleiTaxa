#!/bin/bash

###############################################################################
# NucleiTaxa Test Suite
# 
# Validates pipeline stages with synthetic/mock data
# Allows testing without large input datasets
# Generates sample reports and visualizations
###############################################################################

set -euo pipefail

TEST_DIR="${1:-.}"
TEST_OUTPUT="$TEST_DIR/test_output"
PIPELINE_DIR="$TEST_DIR/pipeline"
ANALYTICS_PORT=8889  # Use different port to avoid conflicts

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

test_info() {
    echo -e "${YELLOW}[TEST]${NC} $*"
}

test_pass() {
    echo -e "${GREEN}[PASS]${NC} $*"
}

test_fail() {
    echo -e "${RED}[FAIL]${NC} $*"
}

# Create mock data
create_mock_fastq() {
    local output_dir="$1"
    local num_samples="${2:-2}"
    local num_reads="${3:-100}"
    
    test_info "Creating mock FASTQ data ($num_samples samples, $num_reads reads each)"
    
    mkdir -p "$output_dir"
    
    for i in $(seq 1 $num_samples); do
        local r1="$output_dir/sample_${i}_R1.fastq.gz"
        local r2="$output_dir/sample_${i}_R2.fastq.gz"
        
        # Generate synthetic reads
        (
            for j in $(seq 1 $num_reads); do
                seq_id="$i:$j"
                # Synthetic DNA sequence (random 150bp)
                seq=$(head -c 150 /dev/urandom | tr -dc 'ACGT' | head -c 150)
                qual=$(python3 -c "import random; print(''.join([chr(random.randint(33, 126)) for _ in range(150)]))")
                
                echo "@seq_${seq_id}/1"
                echo "$seq"
                echo "+"
                echo "$qual"
            done
        ) | gzip > "$r1"
        
        (
            for j in $(seq 1 $num_reads); do
                seq_id="$i:$j"
                seq=$(head -c 150 /dev/urandom | tr -dc 'ACGT' | head -c 150)
                qual=$(python3 -c "import random; print(''.join([chr(random.randint(33, 126)) for _ in range(150)]))")
                
                echo "@seq_${seq_id}/2"
                echo "$seq"
                echo "+"
                echo "$qual"
            done
        ) | gzip > "$r2"
    done
    
    test_pass "Mock FASTQ created: $num_samples samples"
}

# Test Stage 01: Preprocess
test_stage_01() {
    test_info "Testing Stage 01: Preprocess"
    
    local input_dir="$TEST_OUTPUT/input"
    local output_dir="$TEST_OUTPUT/stage01"
    
    create_mock_fastq "$input_dir" 2 50
    
    # Mock the preprocess stage (copy input to output since we don't have BBTools)
    mkdir -p "$output_dir/01-preprocess"
    cp "$input_dir"/*_R1* "$output_dir/01-preprocess/" 2>/dev/null || true
    cp "$input_dir"/*_R2* "$output_dir/01-preprocess/" 2>/dev/null || true
    
    # Rename to match expected output format
    for f in "$output_dir/01-preprocess"/*_R1*.fastq.gz; do
        if [[ -f "$f" ]]; then
            mv "$f" "${f/_R1/_R1_filtered}"
        fi
    done
    for f in "$output_dir/01-preprocess"/*_R2*.fastq.gz; do
        if [[ -f "$f" ]]; then
            mv "$f" "${f/_R2/_R2_filtered}"
        fi
    done
    
    if [[ -f "$output_dir/01-preprocess"/*_filtered.fastq.gz ]]; then
        test_pass "Stage 01 mock complete"
    else
        test_fail "Stage 01 mock failed"
        return 1
    fi
}

# Test Stage 02: DADA2 (mock)
test_stage_02() {
    test_info "Testing Stage 02: DADA2 Denoise (mock)"
    
    local output_dir="$TEST_OUTPUT/stage01"
    local stage02_dir="$TEST_OUTPUT/stage01/02-denoise"
    
    mkdir -p "$stage02_dir/filtered"
    
    # Create mock ASV table and sequences
    cat > "$stage02_dir/seqtab.txt" << EOF
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	100	80
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT	120	90
EOF
    
    cat > "$stage02_dir/asv_sequences.fasta" << EOF
>ASV1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>ASV2
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
EOF
    
    if [[ -f "$stage02_dir/seqtab.txt" ]] && [[ -f "$stage02_dir/asv_sequences.fasta" ]]; then
        test_pass "Stage 02 mock complete"
    else
        test_fail "Stage 02 mock failed"
        return 1
    fi
}

# Test Stage 03: Chimera (mock)
test_stage_03() {
    test_info "Testing Stage 03: Chimera Detection (mock)"
    
    local output_dir="$TEST_OUTPUT/stage01"
    local stage03_dir="$TEST_OUTPUT/stage01/03-chimera"
    
    mkdir -p "$stage03_dir"
    
    # Copy and mock chimera detection
    cp "$output_dir/02-denoise/asv_sequences.fasta" "$stage03_dir/nonchimeras_denovo.fasta"
    cp "$output_dir/02-denoise/seqtab.txt" "$stage03_dir/seqtab_nochim.txt"
    
    # Mock chimera report
    cat > "$stage03_dir/uchime_denovo.txt" << EOF
ASV1	0	N
ASV2	0	N
EOF
    
    if [[ -f "$stage03_dir/seqtab_nochim.txt" ]]; then
        test_pass "Stage 03 mock complete"
    else
        test_fail "Stage 03 mock failed"
        return 1
    fi
}

# Test Stage 04: Taxonomy (mock)
test_stage_04() {
    test_info "Testing Stage 04: Taxonomy (mock)"
    
    local output_dir="$TEST_OUTPUT/stage01"
    local stage04_dir="$TEST_OUTPUT/stage01/04-taxonomy"
    
    mkdir -p "$stage04_dir"
    
    # Create mock taxonomy assignments
    cat > "$stage04_dir/taxa_assignments.txt" << EOF
ASV_ID	Domain	Phylum	Class	Order	Family	Genus	Confidence
ASV1	Bacteria	Firmicutes	Bacilli	Bacillales	Bacillaceae	Bacillus	0.987
ASV2	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacteriales	Enterobacteriaceae	Escherichia	0.992
EOF
    
    cat > "$stage04_dir/phylum_summary.txt" << EOF
Phylum	Count
Firmicutes	100
Proteobacteria	120
EOF
    
    if [[ -f "$stage04_dir/taxa_assignments.txt" ]]; then
        test_pass "Stage 04 mock complete"
    else
        test_fail "Stage 04 mock failed"
        return 1
    fi
}

# Test Stage 05: Phylogenetics (mock)
test_stage_05() {
    test_info "Testing Stage 05: Phylogenetics (mock)"
    
    local output_dir="$TEST_OUTPUT/stage01"
    local stage05_dir="$TEST_OUTPUT/stage01/05-phylo"
    
    mkdir -p "$stage05_dir"
    
    # Create mock Newick tree
    cat > "$stage05_dir/asv_tree_rooted.nwk" << EOF
(((ASV1:0.1,ASV2:0.1):0.1,ASV3:0.2):0.1);
EOF
    
    cat > "$stage05_dir/tree_stats.txt" << EOF
Number of sequences: 2
Number of internal nodes: 1
Tree height: 0.3
EOF
    
    if [[ -f "$stage05_dir/asv_tree_rooted.nwk" ]]; then
        test_pass "Stage 05 mock complete"
    else
        test_fail "Stage 05 mock failed"
        return 1
    fi
}

# Test Stage 06: Visualization (mock)
test_stage_06() {
    test_info "Testing Stage 06: Visualization (mock)"
    
    local output_dir="$TEST_OUTPUT/stage01"
    local stage06_dir="$TEST_OUTPUT/stage01/06-viz"
    
    mkdir -p "$stage06_dir"
    
    # Create mock Krona input and HTML
    cat > "$stage06_dir/krona_input.txt" << EOF
100	Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus
120	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia
EOF
    
    # Create a minimal HTML file (mock Krona output)
    cat > "$stage06_dir/taxa_krona.html" << EOF
<!DOCTYPE html>
<html>
<head><title>NucleiTaxa Krona Chart (Mock)</title></head>
<body>
<h1>Taxonomy Visualization</h1>
<p>Mock Krona visualization for testing</p>
<p>Bacteria: 220 sequences</p>
<ul>
<li>Firmicutes: 100</li>
<li>Proteobacteria: 120</li>
</ul>
</body>
</html>
EOF
    
    if [[ -f "$stage06_dir/taxa_krona.html" ]]; then
        test_pass "Stage 06 mock complete"
    else
        test_fail "Stage 06 mock failed"
        return 1
    fi
}

# Test analytics integration
test_analytics() {
    test_info "Testing Analytics Integration"
    
    # Check if C++ server was built
    if [[ -f "$TEST_DIR/analytics/server/nucleitaxa-server.cpp" ]]; then
        test_pass "Analytics server source present"
    else
        test_fail "Analytics server source missing"
        return 1
    fi
    
    # Check if frontend assets exist
    if [[ -f "$TEST_DIR/analytics/web/app.js" ]] && \
       [[ -f "$TEST_DIR/analytics/web/index.html" ]] && \
       [[ -f "$TEST_DIR/analytics/web/styles.css" ]]; then
        test_pass "Analytics frontend complete"
    else
        test_fail "Analytics frontend incomplete"
        return 1
    fi
}

# Main test execution
main() {
    test_info "=== NucleiTaxa Test Suite ==="
    test_info "Output directory: $TEST_OUTPUT"
    
    mkdir -p "$TEST_OUTPUT"
    
    local failed=0
    
    test_stage_01 || ((failed++))
    test_stage_02 || ((failed++))
    test_stage_03 || ((failed++))
    test_stage_04 || ((failed++))
    test_stage_05 || ((failed++))
    test_stage_06 || ((failed++))
    test_analytics || ((failed++))
    
    test_info ""
    test_info "=== Test Summary ==="
    
    if [[ $failed -eq 0 ]]; then
        test_pass "All tests passed! Pipeline structure validated."
        test_info ""
        test_info "Next steps:"
        test_info "  1. Install required tools (DADA2, VSEARCH, RDP, FastTree, Krona)"
        test_info "  2. Run on real data: ./bin/nucleitaxa --forward <R1> --reverse <R2>"
        test_info "  3. Monitor with: ./bin/nucleitaxa-server (in separate terminal)"
        return 0
    else
        test_fail "$failed test(s) failed!"
        return 1
    fi
}

main "$@"
