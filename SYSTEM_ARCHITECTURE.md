# NucleiTaxa Bash Refactor - Complete System Architecture

## System Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     NucleiTaxa Bash Pipeline                     │
│                 (2025 Cutting-Edge Bioinformatics)               │
└─────────────────────────────────────────────────────────────────┘

┌─── INPUT ───────────────────────────────────────────────────────┐
│ Paired-end FASTQ files (raw sequencing data)                    │
│ - 16S rRNA (default): Illumina MiSeq, NovaSeq, etc.             │
│ - ITS (fungi): Customizable read length                          │
│ - Custom amplicons: Via config file                              │
└─────────────────────────────────────────────────────────────────┘

                              │
                              ▼

┌─── ORCHESTRATION ───────────────────────────────────────────────┐
│                                                                   │
│  bin/nucleitaxa (Bash CLI)                                       │
│  ├─ Argument parsing (--forward, --reverse, --profile, etc.)     │
│  ├─ Stage execution loop                                         │
│  ├─ Error handling & logging                                     │
│  ├─ Resume capability (--resume-from)                            │
│  ├─ Dry-run mode (--dry-run)                                     │
│  └─ Analytics integration (WebSocket client)                     │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘

                    │         │         │         │         │
                    ▼         ▼         ▼         ▼         ▼

┌──────────────┬──────────────┬──────────────┬──────────────┬──┐
│  Stage 01    │  Stage 02    │  Stage 03    │  Stage 04    │..│
│  Preprocess  │  Denoise     │  Chimera     │  Taxonomy    │  │
│              │              │              │              │  │
│ BBTools/     │ DADA2        │ VSEARCH      │ RDP          │  │
│ Cutadapt     │ (R)          │ UCHIME       │ Classifier   │  │
│              │              │              │ (Java)       │  │
│ QC, trim,    │ Error model, │ De novo +    │ Bayesian     │  │
│ filter       │ ASV infer,   │ reference    │ assignment   │  │
│              │ merge,       │ detection    │ w/ bootstrap │  │
│              │ chimera rm   │ (hybrid)     │ confidence   │  │
└──────────────┴──────────────┴──────────────┴──────────────┴──┘

                    │         │         │         │
                    └─────────┼─────────┼─────────┘
                              ▼
                    
                ┌─────────────────────────────┐
                │  Stage 05: Phylogenetics    │
                │  FastTree 2 (ML, GTR+gamma) │
                │  - Unaligned input          │
                │  - 1000x faster than UPGMA  │
                │  - Rooted tree output       │
                └─────────────────────────────┘

                              ▼

                ┌─────────────────────────────┐
                │  Stage 06: Visualization    │
                │  Krona (Interactive HTML)   │
                │  - Hierarchical pie charts  │
                │  - Zooming & exploration    │
                │  - Standalone output        │
                └─────────────────────────────┘

                              │
                              ▼

┌─── ANALYTICS (LIVE) ────────────────────────────────────────────┐
│                                                                   │
│  WebSocket Architecture:                                          │
│  ┌──────────────────────────┬──────────────────────────┐         │
│  │  Backend (C++)           │  Frontend (JavaScript)   │         │
│  ├──────────────────────────┼──────────────────────────┤         │
│  │ nucleitaxa-server.cpp    │ app.js (WebSocket client)│         │
│  │ - Polls stage outputs    │ - Connects to :8888      │         │
│  │ - 500ms intervals        │ - Updates metrics        │         │
│  │ - JSON broadcast to all  │ - D3.js taxonomy pie     │         │
│  │   connected clients      │ - Plotly quality charts  │         │
│  │                          │ - Stage timeline UI      │         │
│  │ Metrics streamed:        │ - Live log console       │         │
│  │ - ASV count              │ - Dark theme responsive  │         │
│  │ - Chimera count          │ - Hover & click interact │         │
│  │ - Mean quality           │                          │         │
│  │ - GC%                    │ Open: http://localhost:  │         │
│  │ - Taxonomy distribution  │         8888             │         │
│  │ - Active stage           │                          │         │
│  │ - Elapsed time           │                          │         │
│  └──────────────────────────┴──────────────────────────┘         │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘

                              │
                              ▼

┌─── OUTPUT ──────────────────────────────────────────────────────┐
│                                                                   │
│ results/                                                          │
│ ├── 01-preprocess/           QC-filtered FASTQ + metrics        │
│ ├── 02-denoise/              ASV table + representative seqs    │
│ ├── 03-chimera/              Chimera-cleaned table              │
│ ├── 04-taxonomy/             Taxonomy assignments per ASV       │
│ ├── 05-phylo/                Rooted phylogenetic tree           │
│ └── 06-viz/                  Interactive Krona HTML             │
│     └── taxa_krona.html      ◄─ OPEN IN BROWSER                │
│                                                                   │
│ Additional files:                                                │
│ - *_REPORT.txt files in each stage directory                    │
│ - logs (*.log) for debugging                                    │
│ - JSON outputs for downstream tools (QIIME2, PhyloSeq)         │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘

┌─── DOWNSTREAM INTEGRATION ──────────────────────────────────────┐
│                                                                   │
│ QIIME2:     qiime tools import from seqtab + taxa_assignments  │
│ PhyloSeq:   Load tree + OTU table + taxonomy in R              │
│ Custom:     All outputs in standard formats (FASTA, TSV, JSON) │
│ R/Biocond:  phyloseq, phylodiv, vegan packages                │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

## Technology Stack

### Core Bioinformatics Tools
- **DADA2** (R package) - ASV inference + error correction
- **VSEARCH** (C binary) - UCHIME chimera detection
- **RDP Classifier** (Java) - Bayesian taxonomy assignment
- **FastTree 2** (C binary) - ML phylogenetic tree inference
- **Krona** (Perl scripts) - Interactive taxonomy visualization
- **BBTools** (Java) - QC and adapter trimming
- **Cutadapt** (Python binary) - Precise primer removal

### Orchestration & Analytics
- **Bash** (shell script) - Pipeline orchestration
- **C++** (compiled) - WebSocket analytics server
- **Vanilla JavaScript** - Dashboard frontend (no framework)
- **D3.js v7** - Interactive visualization
- **Plotly.js** - Quality metrics charts
- **HTML5 WebSocket** - Real-time communication

### Infrastructure (Planned)
- **Docker/Podman** - Container orchestration
- **Docker Compose** - Multi-service coordination
- **Ubuntu 24.04** - Base image

## Data Flow Example

```
Input: sample_R1.fastq.gz (5M reads, 150bp paired-end)
       sample_R2.fastq.gz

Stage 01 (5 min):
  BBTools BBDuk → 4.8M reads pass QC (96% retention)
  Cutadapt      → 4.7M reads after primer removal
  Quality trim  → 4.5M reads >= 50bp
  
Stage 02 (DADA2, 8 min):
  Learn error model (both directions)
  Dereplication  → 2.3M unique sequences
  Sample inference (Bayes)
  Paired-read merge (12bp overlap)
  Consensus chimera removal
  Result: 847 ASVs

Stage 03 (Chimera VSEARCH, 3 sec):
  De novo UCHIME → 21 chimeras detected
  Reference UCHIME → 5 additional chimeras
  Final ASVs: 821 (97% retention)

Stage 04 (RDP Classifier, 10 sec):
  Taxonomic assignment (0.5+ confidence)
  821 ASVs → 34 phyla, 89 classes, 156 orders
  100% assigned to domain (Bacteria/Archaea)

Stage 05 (FastTree, 0.5 sec):
  Multiple sequence alignment (if needed)
  ML tree inference (GTR+gamma)
  Midpoint rooting
  821-leaf tree generated

Stage 06 (Krona, 0.3 sec):
  Convert to Krona format
  Generate HTML visualization
  Interactive pie chart generated

Total: 821 ASVs, phylogenetic tree, interactive chart
Time: ~13 minutes
Memory: 4GB peak (DADA2 learning phase)
```

## Quality Control Gates

Each stage validates its output:

```
Stage 01 → Check: *.fastq.gz files exist, gzip valid
Stage 02 → Check: seqtab.txt has >0 rows, FASTA valid
Stage 03 → Check: seqtab_nochim.txt has fewer rows than input
Stage 04 → Check: taxa_assignments.txt has matching ASV count
Stage 05 → Check: asv_tree_rooted.nwk is valid Newick
Stage 06 → Check: taxa_krona.html >1KB and contains tree data
```

Failed validation stops pipeline with clear error message.

## Configuration Management

### Profiles (Preset Configurations)
- **16s** (default): V3-V4 region, Illumina standard params
- **its**: ITS1/ITS2 fungal, longer region, higher sequence diversity
- **custom**: User-provided config file

### Config File Format
```bash
# Preprocess
TRUNCATE_LENGTH_F=240
TRUNCATE_LENGTH_R=160
MIN_QUALITY=20

# DADA2
MAX_EE_F=2
MAX_EE_R=2
MIN_OVERLAP=12

# Chimera
UCHIME_THRESHOLD=0.85

# Taxonomy
RDP_CONFIDENCE_THRESHOLD=0.5
GENE_TYPE="16srrna"

# Resources
JOBS=8
RDP_MAX_MEMORY=4096
```

## Parallelization Strategy

- **Stage 01**: BBTools respects --threads flag (configurable)
- **Stage 02**: DADA2 multithread learning (R parallel package)
- **Stage 03**: VSEARCH --threads flag
- **Stage 04**: Sequential (RDP doesn't parallelize well)
- **Stage 05**: FastTree respects thread count
- **Stage 06**: Single-threaded (negligible time)

**Bottleneck:** Stage 02 DADA2 error learning (single-threaded in R)
**Recommendation:** Use 16+ cores for full parallelism on other stages

## Error Recovery

### Resume Capability
```bash
# If pipeline interrupted at stage 04:
./bin/nucleitaxa \
    --resume-from 04 \
    --output results

# Resumes from stage 04, skips 01-03
```

### Logging & Debugging
- Each stage logs to `<stage>/$(basename stage).log`
- Main orchestrator logs to `$(output_dir)/nucleitaxa.log`
- Use `tail -f` to monitor real-time
- `--dry-run` shows all commands without executing

## Containerization (Planned)

```dockerfile
FROM ubuntu:24.04

# Install bioinformatics tools
RUN apt update && apt install -y \
    r-base r-cran-ape \
    vsearch fasttree kronatools \
    git curl

# DADA2, RDP, etc.
RUN ...

COPY bin/nucleitaxa /usr/local/bin/
COPY pipeline /opt/nucleitaxa/pipeline
COPY analytics /opt/nucleitaxa/analytics

ENTRYPOINT ["/usr/local/bin/nucleitaxa"]
```

## Performance Optimization Tips

1. **Use SSD for temp files** – I/O bottleneck in preprocessing
2. **Allocate 4GB+ RAM** – DADA2 memory requirements
3. **Max CPU cores** – Use `--jobs $(nproc)` for full parallelism
4. **Monitor with analytics** – Real-time insights into bottlenecks
5. **Batch processing** – Submit multiple samples to run in parallel

## Future Extensions

- **SLURM integration** – Job submission to HPC clusters
- **Batch framework** – Automatic multi-sample coordination
- **Database auto-fetch** – RDP, SILVA, Greengenes auto-download
- **Extended output formats** – Excel reports, PDF figures
- **Web UI** – Remote dashboard access
- **Cloud deployment** – AWS/GCP ready containers

---

**System designed for transparency, reproducibility, and 2025 bioinformatics best practices.**
