"""
Quality Control Module

Implements DADA2-inspired quality filtering and denoising.
Pure Python implementation to eliminate R/dada2 dependency.

Bioinformatic principles:
- Quality filtering removes low-confidence bases before clustering
- Expected error calculation weighs uncertainty appropriately
- Chimera detection identifies PCR artifacts with high sensitivity
"""

from typing import List, Tuple, Dict
import numpy as np
from .io import FastqRecord, FastqReader, FastqFilter


class QualityController:
    """
    Orchestrates quality control pipeline for amplicon sequences.
    
    Pipeline steps:
    1. Quality profile analysis (pre-filtering diagnostics)
    2. Filter by quality thresholds (EE-based)
    3. Trim primers and adapters
    4. Dereplicate (group identical sequences)
    5. Optional: Remove chimeras
    """
    
    def __init__(self, config_profile):
        """
        Args:
            config_profile: NuclidProfile instance with QC parameters
        """
        self.profile = config_profile
        self.filter = FastqFilter(
            trim_left=config_profile.trim_left,
            trim_right=config_profile.trim_right,
            trunc_len=(config_profile.trunc_len_f, config_profile.trunc_len_r),
            max_ee=config_profile.max_ee,
            min_quality=config_profile.min_quality,
            max_n=config_profile.max_n,
        )
        self.stats = {
            'input': 0,
            'passed_qc': 0,
            'failed_length': 0,
            'failed_ee': 0,
            'failed_n': 0,
        }
    
    def filter_fastq_pair(self, fwd_file: str, rev_file: str,
                         output_prefix: str) -> Tuple[str, str]:
        """
        Filter paired-end FASTQ files.
        
        Args:
            fwd_file: Forward reads (R1)
            rev_file: Reverse reads (R2)
            output_prefix: Output file prefix
        
        Returns: Tuple of (filtered_fwd, filtered_rev) file paths
        """
        fwd_reader = FastqReader(fwd_file)
        rev_reader = FastqReader(rev_file)
        
        output_fwd = f"{output_prefix}_R1_filtered.fastq.gz"
        output_rev = f"{output_prefix}_R2_filtered.fastq.gz"
        
        # Open output files
        from .io import FastaWriter
        import gzip
        
        with gzip.open(output_fwd, 'wt') as fwd_out, \
             gzip.open(output_rev, 'wt') as rev_out:
            
            for fwd_record, rev_record in zip(fwd_reader, rev_reader):
                fwd_filtered = self.filter.filter_record(fwd_record, is_forward=True)
                rev_filtered = self.filter.filter_record(rev_record, is_forward=False)
                
                if fwd_filtered and rev_filtered:
                    # Write FASTQ format
                    self._write_fastq_record(fwd_out, fwd_filtered)
                    self._write_fastq_record(rev_out, rev_filtered)
                    self.stats['passed_qc'] += 1
                else:
                    if not fwd_filtered or not rev_filtered:
                        self.stats['failed_ee'] += 1
                
                self.stats['input'] += 1
        
        return output_fwd, output_rev
    
    @staticmethod
    def _write_fastq_record(handle, record: FastqRecord):
        """Write record in FASTQ format."""
        handle.write(f"@{record.id}\n{record.seq}\n+\n{record.qual}\n")
    
    def dereplicate(self, fastq_file: str) -> Dict[str, int]:
        """
        Collapse identical sequences (dereplication).
        
        Bioinformatic rationale:
        - Reduces redundancy: identical sequences are counted, not processed separately
        - Speeds up downstream analyses (ASV/OTU clustering)
        - Preserves abundance information
        
        Args:
            fastq_file: Input FASTQ file
        
        Returns: Dict mapping sequence -> count
        """
        derep = {}
        reader = FastqReader(fastq_file)
        
        for record in reader:
            seq = record.seq
            derep[seq] = derep.get(seq, 0) + 1
        
        return derep
    
    def detect_chimeras_denovo(self, sequences: List[str],
                               abundance: List[int],
                               threshold: float = 0.85) -> np.ndarray:
        """
        Detect chimeric sequences *de novo* using abundance ratio.
        
        Bioinformatic principle:
        - True variants arise from biological diversity
        - Chimeras are PCR artifacts: combinations of two parent sequences
        - Chimeric sequences typically have lower abundance than parents
        - Two-way sequence dissimilarity can identify chimera points
        
        Simple heuristic: Flag sequences with suspiciously low abundance
        relative to more similar sequences.
        
        Args:
            sequences: List of DNA sequences
            abundance: Abundance counts for each sequence
            threshold: Minimum sequence similarity to "parent" sequence
        
        Returns: Boolean array (True = chimera, False = non-chimera)
        """
        n = len(sequences)
        is_chimera = np.zeros(n, dtype=bool)
        
        for i, seq_i in enumerate(sequences):
            # Find most abundant sequence with >threshold similarity
            best_similarity = 0
            best_j = -1
            
            for j, seq_j in enumerate(sequences):
                if i == j:
                    continue
                
                sim = self._similarity(seq_i, seq_j)
                if sim > threshold and abundance[j] > abundance[i]:
                    if sim > best_similarity:
                        best_similarity = sim
                        best_j = j
            
            # If found high-similarity parent with higher abundance, likely chimera
            if best_j >= 0:
                is_chimera[i] = True
        
        return is_chimera
    
    @staticmethod
    def _similarity(seq1: str, seq2: str) -> float:
        """Compute sequence similarity (ignoring length differences)."""
        if len(seq1) != len(seq2):
            return 0.0
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)
