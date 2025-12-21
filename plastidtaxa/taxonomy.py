"""
Taxonomic Assignment Module

Classifies sequences against reference databases using k-mer matching,
similarity-based, or machine learning methods.
"""

from typing import List, Dict, Tuple, Optional
import numpy as np
from pathlib import Path


class TaxonomyAssigner:
    """
    Assign sequences to taxonomic lineages using reference database.
    
    Methods:
    1. Exact match: Compare to database, take top hits
    2. k-mer based: Count shared k-mers (fast, scalable)
    3. Machine learning: Train classifier on known sequences
    
    Output: (sequence, kingdom, phylum, class, order, family, genus, species)
    """
    
    def __init__(self, reference_db: str, method: str = "kmer"):
        """
        Args:
            reference_db: FASTA file with reference sequences & taxonomy
            method: "exact", "kmer", or "blast"
        """
        self.reference_db = Path(reference_db)
        self.method = method
        self.references = self._load_references()
        self.kmer_size = 8
    
    def _load_references(self) -> Dict[str, Tuple[str, str]]:
        """
        Load reference sequences and taxonomy.
        
        Format: 
            >seq_id;tax=k__Kingdom;p__Phylum;c__Class;...
            ACGTACGTACGT...
        
        Returns: Dict of {sequence -> (seq_id, taxonomy_string)}
        """
        references = {}
        
        try:
            from Bio import SeqIO
            for record in SeqIO.parse(self.reference_db, "fasta"):
                tax = self._extract_taxonomy(record.description)
                references[str(record.seq).upper()] = (record.id, tax)
        except ImportError:
            print("Warning: Biopython not available, using simple FASTA parser")
            # Fallback: simple line-based parser
            with open(self.reference_db) as f:
                current_seq = None
                current_id = None
                for line in f:
                    if line.startswith('>'):
                        if current_seq:
                            current_id_upper = current_id.upper()
                            references[current_seq] = (current_id_upper, current_id)
                        current_id = line[1:].strip()
                        current_seq = ""
                    else:
                        current_seq += line.strip()
        
        print(f"Loaded {len(references)} reference sequences")
        return references
    
    @staticmethod
    def _extract_taxonomy(description: str) -> str:
        """Extract taxonomy string from FASTA header."""
        # Handle SILVA format: >GBXX123.1 Bacteria;Firmicutes;...
        if ';' in description:
            return description.split(' ', 1)[-1]
        # Handle other formats
        return description
    
    def assign(self, sequence: str, confidence_threshold: float = 0.8) -> Dict:
        """
        Assign taxonomy to a sequence.
        
        Args:
            sequence: DNA sequence (will be converted to uppercase)
            confidence_threshold: Min score (0-1) to accept assignment
        
        Returns: Dict with keys:
            - 'sequence_id': input sequence
            - 'taxonomy': assigned lineage
            - 'confidence': match score (0-1)
            - 'reference_id': matched reference sequence
        """
        seq_upper = sequence.upper()
        
        if self.method == "kmer":
            return self._assign_kmer(seq_upper, confidence_threshold)
        elif self.method == "exact":
            return self._assign_exact(seq_upper, confidence_threshold)
        else:
            raise ValueError(f"Unknown method: {self.method}")
    
    def _assign_kmer(self, sequence: str, threshold: float) -> Dict:
        """
        k-mer based taxonomic assignment.
        
        Bioinformatic principle:
        - Extract all k-mers from query sequence
        - Count shared k-mers with each reference
        - Rank references by similarity
        - Return top match if confidence >= threshold
        """
        query_kmers = self._extract_kmers(sequence)
        if not query_kmers:
            return {"taxonomy": "unclassified", "confidence": 0.0}
        
        best_score = 0
        best_ref = None
        best_tax = None
        
        for ref_seq, (ref_id, tax) in self.references.items():
            ref_kmers = self._extract_kmers(ref_seq)
            shared = len(query_kmers & ref_kmers)
            total = len(query_kmers | ref_kmers)
            
            if total == 0:
                continue
            
            jaccard_sim = shared / total
            
            if jaccard_sim > best_score:
                best_score = jaccard_sim
                best_ref = ref_id
                best_tax = tax
        
        confidence = best_score if best_score >= threshold else 0.0
        
        return {
            "taxonomy": best_tax or "unclassified",
            "confidence": confidence,
            "reference_id": best_ref or "unknown",
        }
    
    def _assign_exact(self, sequence: str, threshold: float) -> Dict:
        """Exact match assignment."""
        if sequence in self.references:
            ref_id, tax = self.references[sequence]
            return {
                "taxonomy": tax,
                "confidence": 1.0,
                "reference_id": ref_id,
            }
        
        # Fallback: find closest by Hamming distance
        best_dist = float('inf')
        best_ref = None
        best_tax = None
        
        for ref_seq, (ref_id, tax) in self.references.items():
            if len(ref_seq) == len(sequence):
                dist = sum(1 for a, b in zip(sequence, ref_seq) if a != b)
                if dist < best_dist:
                    best_dist = dist
                    best_ref = ref_id
                    best_tax = tax
        
        if best_dist > float('inf'):
            similarity = max(0, 1.0 - (best_dist / len(sequence)))
            if similarity >= threshold:
                return {
                    "taxonomy": best_tax,
                    "confidence": similarity,
                    "reference_id": best_ref,
                }
        
        return {"taxonomy": "unclassified", "confidence": 0.0}
    
    def _extract_kmers(self, sequence: str, k: Optional[int] = None) -> set:
        """Extract all k-mers from sequence."""
        if k is None:
            k = self.kmer_size
        
        kmers = set()
        for i in range(len(sequence) - k + 1):
            kmers.add(sequence[i:i+k])
        return kmers
    
    def batch_assign(self, sequences: List[str]) -> List[Dict]:
        """Assign taxonomy to multiple sequences."""
        results = []
        for seq in sequences:
            results.append(self.assign(seq))
        return results
