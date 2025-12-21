"""
Denoising & Clustering Module

Pure Python implementations of ASV (Amplicon Sequence Variant) inference.
Replaces DADA2 R package with open-source Python alternatives.
"""

from typing import List, Dict, Tuple
import numpy as np
from collections import Counter
from difflib import SequenceMatcher


class Denoiser:
    """
    Denoise amplicon sequences to identify true variants (ASVs).
    
    Bioinformatic principle:
    - Single-pass clustering: Group reads by edit distance
    - Abundance-based filtering: Low-abundance clusters likely errors
    - Error model: Base substitution errors correlate with quality scores
    - Output: Amplicon Sequence Variants (ASVs) with abundance counts
    """
    
    def __init__(self, error_threshold: float = 0.01):
        """
        Args:
            error_threshold: Max edit distance for clustering (0-1)
        """
        self.error_threshold = error_threshold
    
    def cluster_sequences(self, sequences: List[str],
                         abundance: List[int],
                         min_abundance: int = 2) -> Tuple[List[str], List[int]]:
        """
        Cluster sequences by similarity to identify true variants.
        
        Algorithm: Greedy clustering
        - Sort by decreasing abundance (most abundant = most likely true)
        - Compare each to existing clusters
        - Assign to first similar cluster or create new
        
        Args:
            sequences: List of DNA sequences
            abundance: Abundance count for each sequence
            min_abundance: Minimum reads to retain (filter singletons)
        
        Returns: (clustered_seqs, clustered_abundance)
        """
        # Filter by minimum abundance
        filtered = [(s, a) for s, a in zip(sequences, abundance) if a >= min_abundance]
        if not filtered:
            return [], []
        
        sequences, abundance = zip(*sorted(filtered, key=lambda x: -x[1]))
        
        clusters = {}  # sequence -> total_abundance
        assignments = {}  # sequence_idx -> cluster_representative
        
        for seq, abund in zip(sequences, abundance):
            assigned = False
            
            # Try to match to existing cluster
            for cluster_rep in clusters.keys():
                if self._similarity(seq, cluster_rep) >= (1.0 - self.error_threshold):
                    clusters[cluster_rep] += abund
                    assignments[seq] = cluster_rep
                    assigned = True
                    break
            
            # Create new cluster
            if not assigned:
                clusters[seq] = abund
                assignments[seq] = seq
        
        return list(clusters.keys()), list(clusters.values())
    
    @staticmethod
    def _similarity(seq1: str, seq2: str) -> float:
        """Compute Levenshtein-inspired similarity."""
        if len(seq1) != len(seq2):
            return 0.0
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)
    
    def error_model(self, qual_string: str, seq_position: int) -> float:
        """
        Estimate error probability at a sequence position.
        
        Bioinformatic principle:
        - Quality scores (Phred) inversely correlate with error probability
        - Q = -10 * log10(P_error)
        - Rearranged: P_error = 10^(-Q/10)
        
        Args:
            qual_string: Quality scores (Phred+33 ASCII)
            seq_position: 0-indexed base position
        
        Returns: Estimated error probability (0-1)
        """
        if seq_position >= len(qual_string):
            return 1.0
        
        q = ord(qual_string[seq_position]) - 33
        return 10 ** (-q / 10)


class SequenceMerger:
    """Merge paired-end reads into consensus sequences."""
    
    def __init__(self, min_overlap: int = 10, max_mismatch: int = 0):
        """
        Args:
            min_overlap: Minimum overlap required to merge
            max_mismatch: Max mismatches in overlap region
        """
        self.min_overlap = min_overlap
        self.max_mismatch = max_mismatch
    
    def merge_pair(self, fwd_seq: str, rev_seq: str,
                   fwd_qual: str, rev_qual: str) -> Tuple[str, str, bool]:
        """
        Merge forward and reverse reads into consensus sequence.
        
        Bioinformatic logic:
        - Find overlap between 3' end of forward and 5' end of reverse complement
        - Validate sufficient overlap (min_overlap)
        - For disagreements, choose base from higher quality score
        - Return merged sequence or None if merge fails
        
        Args:
            fwd_seq, rev_seq: DNA sequences
            fwd_qual, rev_qual: Quality strings
        
        Returns: (merged_seq, merged_qual, success)
        """
        from Bio.Seq import reverse_complement
        
        rev_rc = str(reverse_complement(rev_seq))
        rev_qual_rc = rev_qual[::-1]
        
        # Find best overlap
        best_overlap = 0
        best_overlap_pos = -1
        
        for shift in range(1, min(len(fwd_seq), len(rev_rc)) + 1):
            fwd_end = fwd_seq[-shift:]
            rev_start = rev_rc[:shift]
            
            mismatches = sum(1 for a, b in zip(fwd_end, rev_start) if a != b)
            if mismatches <= self.max_mismatch and shift >= self.min_overlap:
                best_overlap = shift
                best_overlap_pos = len(fwd_seq) - shift
        
        if best_overlap == 0:
            return "", "", False
        
        # Build consensus from overlap region
        merged = fwd_seq[:best_overlap_pos]
        merged_qual = fwd_qual[:best_overlap_pos]
        
        for i in range(best_overlap):
            fwd_pos = best_overlap_pos + i
            rev_pos = i
            
            fwd_q = ord(fwd_qual[fwd_pos]) - 33 if fwd_pos < len(fwd_qual) else 0
            rev_q = ord(rev_qual_rc[rev_pos]) - 33 if rev_pos < len(rev_qual_rc) else 0
            
            if fwd_q > rev_q:
                merged += fwd_seq[fwd_pos]
                merged_qual += fwd_qual[fwd_pos]
            else:
                merged += rev_rc[rev_pos]
                merged_qual += rev_qual_rc[rev_pos]
        
        # Add remaining 3' overhang from reverse complement
        remaining_start = best_overlap
        merged += rev_rc[remaining_start:]
        merged_qual += rev_qual_rc[remaining_start:]
        
        return merged, merged_qual, True
