# PlastidTaxa Refactor – Architecture & Implementation Guide

## Overview

This refactored version transforms PlastidTaxa from R-dependent scripts into a **modular, scalable Python package** supporting multiple nucleid targets (plastid 23S/16S, mitochondrial, ITS, 16S, etc.).

### Key Improvements

| Aspect | Before | After |
|--------|--------|-------|
| **Architecture** | Monolithic R scripts | Modular Python library |
| **Nucleid Support** | Hardcoded for plastid 23S | Configurable profiles (5+ targets) |
| **Dependencies** | R/dada2, Python, rpy2 | Pure Python (pip installable) |
| **Scalability** | Fixed pipeline | Pluggable stages & methods |
| **Security** | Biopython ≤1.86 (XXE) | Biopython ≥1.87 + defusedxml |
| **Readability** | Minimal comments | Inline bioinformatic documentation |
| **Installation** | Manual setup | `pip install -e .` |

---

## Architecture: Module Overview

```
plastidtaxa/
├── __init__.py          # Main package exports
├── config.py            # Nucleid profiles (23S, 16S, ITS, etc.)
├── io.py                # FASTQ/FASTA reading, quality metrics
├── fetch.py             # NCBI Entrez (XXE-safe), database downloads
├── qc.py                # Quality control, filtering, dereplication
├── denoise.py           # ASV clustering, pair merging
├── taxonomy.py          # k-mer & exact-match classification
├── phylo.py             # Distance matrices, neighbor-joining trees
├── cli.py               # Command-line interface
```

### Module Dependencies

```
fetch.py (XXE-safe Entrez)
         ↓
config.py (profiles)
    ↓       ↓
   io.py   qc.py (filtering, dereplication, chimera detection)
    ↓       ↓
denoise.py (ASV inference)
    ↓
taxonomy.py (k-mer classification)
    ↓
phylo.py (tree construction)
    ↓
    cli.py (orchestration)
```

---

## Installation & Usage

### Installation

```bash
# Clone and install
git clone https://github.com/xaoscience/PlastidTaxa.git
cd PlastidTaxa
pip install -e .

# Or install with dev tools
pip install -e ".[dev]"

# Optional: Jupyter notebook support
pip install -e ".[notebook]"
```

### CLI Usage

```bash
# List available profiles
plastidtaxa list-profiles

# Run full pipeline (plastid 23S)
plastidtaxa run-pipeline \
    --profile plastid_23s \
    --forward reads_R1.fastq.gz \
    --reverse reads_R2.fastq.gz \
    --reference silva_23s.fasta \
    --output ./results

# Run pipeline (fungal ITS)
plastidtaxa run-pipeline \
    --profile fungal_its \
    --forward reads_R1.fastq.gz \
    --reverse reads_R2.fastq.gz \
    --reference unite_database.fasta \
    --output ./results_ITS

# Generate QC report
plastidtaxa qc-report --input reads_R1.fastq.gz --profile plastid_23s
```

### Python API

```python
from plastidtaxa import get_profile, QualityController, Denoiser, TaxonomyAssigner

# Get profile
profile = get_profile("plastid_23s")

# Quality control
qc = QualityController(profile)
fwd_filt, rev_filt = qc.filter_fastq_pair("R1.fastq.gz", "R2.fastq.gz", "filtered")

# Dereplication
derep = qc.dereplicate("filtered_R1_filtered.fastq.gz")

# Denoising
denoiser = Denoiser()
asv_seqs, asv_counts = denoiser.cluster_sequences(
    list(derep.keys()),
    list(derep.values())
)

# Taxonomy
assigner = TaxonomyAssigner("ref_db.fasta", method="kmer")
results = assigner.batch_assign(asv_seqs)
```

---

## Bioinformatic Principles (Documented in Code)

### 1. **Quality Control** (`qc.py`)

**Expected Error (EE) Filtering:**
- Formula: $EE = \sum_{i=1}^{n} 10^{-Q_i/10}$
- Why: Phred scores weight error probability; sum of errors predicts sequence reliability
- Default threshold: EE ≤ 2.0 (high confidence)

**Chimera Detection:**
- Flag sequences with high similarity to more-abundant sequences
- Heuristic: If `abundance[A] > abundance[B]` and `similarity(A, B) > 0.85`, then B is likely chimeric
- Future: Implement UCHIME or vsearch for more robust detection

### 2. **Denoising** (`denoise.py`)

**ASV Clustering:**
- Group identical/near-identical sequences by edit distance
- Retain abundance information (unlike traditional OTU clustering)
- More sensitive to rare variants than 97% OTU clustering

**Pair Merging:**
- Overlap forward (R1) & reverse complement (R2) sequences
- Consensus base selection: Choose base with higher quality score
- Minimum overlap: 10 bp (configurable)

### 3. **Taxonomic Assignment** (`taxonomy.py`)

**k-mer Method:**
- Extract all k-mers (k=8 by default) from query
- Compute Jaccard similarity: $J = \frac{|A \cap B|}{|A \cup B|}$ vs. reference
- Rank references; report top hit if confidence ≥ threshold

**Exact Match Method:**
- Direct lookup; fallback to Hamming distance if no exact match
- Converts Hamming distance to similarity: $sim = 1 - \frac{H}{L}$

### 4. **Phylogenetics** (`phylo.py`)

**Jukes-Cantor Distance:**
- Formula: $d = -\frac{3}{4} \ln(1 - \frac{4p}{3})$ where $p$ = fraction of differences
- Assumption: All substitutions equally likely
- Use for distant sequences or when mutation model unknown

**Kimura 2-Parameter:**
- Distinguishes transitions (A↔G, C↔T) from transversions
- More realistic; use for closely related sequences
- Formula handles unequal transition/transversion ratios

**Neighbor-Joining:**
- Agglomerative clustering: Repeatedly join closest pairs
- Corrects distances for unequal divergence rates
- More accurate than UPGMA for realistic evolutionary rates

---

## Security: XXE Vulnerability Mitigation

### The Vulnerability

Biopython ≤ 1.86 `Bio.Entrez.read()` uses unsafe XML parsing:

```python
# OLD (vulnerable)
from Bio import Entrez
handle = Entrez.efetch(db="nucleotide", id="123456", rettype="gb", retmode="xml")
record = Entrez.read(handle)  # Vulnerable to XXE injection!
```

Attackers could inject malicious DTD/entities to:
- Exfiltrate local files
- Perform SSRF attacks
- Cause DoS via billion laughs attack

### Our Solution

**Option 1: Upgrade Biopython (Recommended)**
```bash
pip install biopython>=1.87
```

**Option 2: Use defusedxml (Fallback)**
```python
from plastidtaxa.fetch import safe_entrez_fetch

# XXE-safe wrapper
results = safe_entrez_fetch("16S ribosomal RNA[Title]", db="nucleotide")
```

Implementation in [fetch.py](plastidtaxa/fetch.py):
- Uses `defusedxml.ElementTree` instead of standard `xml.etree`
- Disables DTD processing and external entity resolution
- Rate limits NCBI requests (0.5 sec delay between calls)

---

## Extending the Pipeline

### Add a New Nucleid Profile

```python
# In config.py, add to PROFILES dict:
"custom_target": NuclidProfile(
    name="custom_target",
    description="Description of your target",
    region="Custom_Region",
    expected_length=1000,
    trim_left=(10, 10),
    trim_right=(0, 0),
    trunc_len_f=250,
    trunc_len_r=200,
    max_ee=(2.0, 2.0),
    min_quality=3,
    max_n=0,
)
```

Then use:
```bash
plastidtaxa run-pipeline --profile custom_target ...
```

### Implement a New Taxonomy Method

```python
# In taxonomy.py, extend TaxonomyAssigner:
def _assign_custom(self, sequence: str, threshold: float) -> Dict:
    """Your custom assignment logic."""
    # Implement...
    return {"taxonomy": "...", "confidence": ...}

# Update __init__ to add method option
```

### Add Advanced Filtering

```python
# In qc.py, extend FastqFilter:
def filter_by_gc_content(self, seq: str, gc_min=0.4, gc_max=0.6) -> bool:
    """Filter sequences by GC content."""
    gc = (seq.count('G') + seq.count('C')) / len(seq)
    return gc_min <= gc <= gc_max
```

---

## Comparison: R DADA2 vs. Python Implementation

| Feature | DADA2 (R) | PlastidTaxa (Python) |
|---------|-----------|----------------------|
| **Speed** | Fast (compiled C/C++) | Slower, pure Python |
| **Accuracy** | Advanced error model | Simplified heuristic |
| **Parallelization** | Built-in (`multithread=TRUE`) | Easily parallelizable (numpy/dask) |
| **Chimera Detection** | UCHIME algorithm | Simple abundance-based |
| **Installation** | Complex (R + bioconda) | Simple (`pip install`) |
| **Portability** | Works on all platforms | Python 3.10+ requirement |

**Recommendation:**
- Use `dada2` R/Python package if maximum sensitivity required
- Use PlastidTaxa for quick, containerized, scalable pipelines

---

## Testing & Validation

```bash
# Run tests (placeholder - add pytest tests)
pytest tests/

# Run example pipeline
plastidtaxa run-pipeline \
    --profile plastid_23s \
    --forward data/example_R1.fastq.gz \
    --reverse data/example_R2.fastq.gz \
    --reference data/silva_23s_subset.fasta \
    --output ./test_output
```

---

## Future Enhancements

- [ ] Integrate UCHIME for robust chimera detection
- [ ] Add Bayesian classifier option (similar to RDP)
- [ ] Implement VSEARCH wrapper for clustering
- [ ] Parallel processing (multiprocessing/dask)
- [ ] Visualization module (matplotlib/plotly)
- [ ] Support for single-end reads
- [ ] Database auto-download (SILVA, UNITE, RDP)
- [ ] Jupyter notebook workflows

---

## References

1. **DADA2:** Callahan et al. (2016). "Exact sequence variants ecologically reflect bacterial community structure." *mSystems* 1(3).
2. **XXE Prevention:** OWASP. "XML External Entity (XXE) Prevention Cheat Sheet."
3. **Phylogenetics:** Saitou & Nei (1987). "The neighbor-joining method: a new method for reconstructing phylogenetic trees."
4. **Defusedxml:** Viraptor. "The dangers of XML parsing."

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

GPL-3.0. See [LICENSE](LICENSE).
