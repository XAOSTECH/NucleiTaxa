"""
PlastidTaxa Core Module

Nucleid-agnostic taxonomy pipeline with modular architecture.
Supports plastid (23S, 16S), mitochondrial, fungal ITS, and other targets.

Architecture:
  - io.py: Data I/O (FASTQ/FASTA/GenBank)
  - fetch.py: Remote sequence retrieval (NCBI Entrez, XXE-safe)
  - qc.py: Quality control & filtering
  - denoise.py: Denoising & ASV inference
  - taxonomy.py: Taxonomic assignment
  - phylo.py: Phylogenetic tree construction
  - config.py: Configuration & nucleid profiles
"""

__version__ = "0.2.0"
__author__ = "xaoscience"

from .config import NuclidProfile, get_profile
from .io import FastqReader, FastaWriter
from .fetch import safe_entrez_fetch
from .qc import QualityController
from .denoise import Denoiser
from .taxonomy import TaxonomyAssigner
from .phylo import PhyloBuilder

__all__ = [
    "NuclidProfile",
    "get_profile",
    "FastqReader",
    "FastaWriter",
    "safe_entrez_fetch",
    "QualityController",
    "Denoiser",
    "TaxonomyAssigner",
    "PhyloBuilder",
]
