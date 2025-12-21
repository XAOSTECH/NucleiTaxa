# REFACTOR BRANCH STATUS

## Summary

Created comprehensive refactor of PlastidTaxa from monolithic R scripts to modular, scalable Python package.

### What Changed

✅ **Package Structure**
- Converted to pip-installable `setup.py` + `pyproject.toml`-ready structure
- Added `plastidtaxa/` module with 8 focused files (io, config, qc, denoise, taxonomy, phylo, fetch, cli)

✅ **Scalability**
- Nucleid-agnostic design: 5+ pre-configured profiles (plastid_23s, plastid_16s, mito_12s, fungal_its, bacterial_16s)
- Pluggable methods for taxonomy assignment (k-mer, exact match, extensible)
- Profile system in `config.py` allows users to add custom targets

✅ **Security: XXE Vulnerability Fixed**
- Replaced Biopython ≤1.86 with ≥1.87 (patched)
- Added `defusedxml` for safe XML parsing in `fetch.py`
- Implemented `safe_entrez_fetch()` wrapper preventing XXE injection

✅ **Dependency Reduction**
- **Removed:** R, rpy2, dada2 (R package)
- **Replaced with:** Pure Python implementations
  - QC/filtering: Custom expected-error algorithm
  - Denoising: Greedy clustering for ASV inference
  - Phylo: Neighbor-joining tree construction
- Minimal Python deps: numpy, scipy, pandas, biopython, defusedxml, click

✅ **Modularity & Documentation**
- Each module has clear responsibility (io, config, qc, denoise, taxonomy, phylo, fetch)
- Inline bioinformatic explanations: Why algorithms work, formulas, assumptions
- Examples in docstrings for all major classes

✅ **CLI Interface**
- New `plastidtaxa` command-line tool (via Click framework)
- Commands: `list-profiles`, `run-pipeline`, `qc-report`
- Orchestrates full pipeline: QC → dereplication → denoising → taxonomy → output

✅ **Documentation**
- `REFACTOR_GUIDE.md`: Comprehensive architecture guide
- Installation, usage (CLI + Python API), bioinformatic principles
- Security discussion (XXE), testing, extension points

### Files Added/Modified

**New:**
- `setup.py` – Package configuration
- `plastidtaxa/__init__.py` – Module exports
- `plastidtaxa/config.py` – Nucleid profiles
- `plastidtaxa/io.py` – FASTQ/FASTA I/O + quality metrics
- `plastidtaxa/fetch.py` – XXE-safe Entrez fetching
- `plastidtaxa/qc.py` – Quality control & filtering
- `plastidtaxa/denoise.py` – ASV clustering & merging
- `plastidtaxa/taxonomy.py` – k-mer & exact-match classification
- `plastidtaxa/phylo.py` – Phylogenetic tree construction
- `plastidtaxa/cli.py` – Command-line interface
- `REFACTOR_GUIDE.md` – Implementation guide

**Updated:**
- `requirements.txt` – Remove R, add defusedxml, bump Biopython
- (Legacy R scripts remain in `scripts/` for reference)

### XXE Fix Details

**Before:** 
```python
from Bio import Entrez
record = Entrez.read(handle)  # Vulnerable in ≤1.86
```

**After:**
```python
from plastidtaxa.fetch import safe_entrez_fetch
results = safe_entrez_fetch("query")  # XXE-protected
```

- Uses `defusedxml.ElementTree` instead of unsafe `xml.etree`
- Disables DTD processing
- Rate limits NCBI (0.5s delay per request)

### Testing Recommendations

1. **Unit tests** (add to `tests/` directory):
   - Config profile loading
   - I/O (FASTQ reading, quality calculation)
   - QC filtering (expected error, N counts)
   - Denoising (clustering, merging)
   - Taxonomy (k-mer matching)
   - Phylo (distance matrix, NJ tree)

2. **Integration test**:
   ```bash
   plastidtaxa run-pipeline \
       --profile plastid_23s \
       --forward data/test_R1.fastq.gz \
       --reverse data/test_R2.fastq.gz \
       --reference data/silva_subset.fasta \
       --output test_results/
   ```

3. **Validation**:
   - Compare output tables to legacy R pipeline on subset data
   - Verify ASV counts, taxonomy assignments, tree topology

### Next Steps

1. **Test on real data** – Validate results match (or exceed) legacy pipeline
2. **Performance optimization** – Profile for bottlenecks, consider parallelization (dask/multiprocessing)
3. **Implement missing features**:
   - UCHIME chimera detection (currently basic abundance-based)
   - VSEARCH/BLAST integration
   - Bayesian/RDP-style classifier
4. **Add pytest test suite** – Cover all modules
5. **Jupyter notebook examples** – Showcase API usage
6. **Docker/Singularity** – Container for reproducibility

### Backward Compatibility

- **Not fully compatible** – API is different (Python vs. R)
- **Legacy scripts** remain in `scripts/` for reference
- Users should use new Python API or CLI

### Branch: `refactor`

To merge to `main` when ready:
```bash
git checkout main
git merge refactor --no-ff -m "Refactor: Modular Python architecture with XXE fix & scalability"
```

---

**Status:** ✅ **Ready for testing & validation**
