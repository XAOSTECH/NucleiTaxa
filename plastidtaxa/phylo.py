"""
Phylogenetic Tree Construction Module

Build neighbor-joining or UPGMA trees from aligned sequences.
Replaces DECIPHER R package with pure Python methods.
"""

from typing import List, Dict, Tuple
import numpy as np
from pathlib import Path


class PhyloBuilder:
    """
    Construct phylogenetic trees from aligned DNA sequences.
    
    Methods:
    1. Distance matrix: Compute pairwise sequence distances
    2. Tree inference: Neighbor-joining or UPGMA clustering
    3. Tree I/O: Read/write Newick format
    
    Bioinformatic principle:
    - Phylogenetic reconstruction groups similar sequences (evolution/descent)
    - Distance metrics account for evolutionary model (Jukes-Cantor, etc.)
    - Tree topology reflects evolutionary relationships
    """
    
    def __init__(self, model: str = "jukes_cantor"):
        """
        Args:
            model: Distance model ("jukes_cantor", "kimura2p")
        """
        self.model = model
    
    def distance_matrix(self, sequences: List[str]) -> np.ndarray:
        """
        Compute pairwise distance matrix from aligned sequences.
        
        Args:
            sequences: List of aligned DNA sequences (same length)
        
        Returns: (n x n) distance matrix
        """
        n = len(sequences)
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                dist = self._compute_distance(sequences[i], sequences[j])
                distances[i, j] = dist
                distances[j, i] = dist
        
        return distances
    
    def _compute_distance(self, seq1: str, seq2: str) -> float:
        """
        Compute evolutionary distance between two sequences.
        
        Jukes-Cantor model:
        - Assumes all substitutions equally likely
        - d = -0.75 * ln(1 - 4p/3)  where p = fraction of differences
        
        Kimura 2-parameter:
        - Distinguishes transitions (A<->G, C<->T) from transversions
        - More realistic for molecular data
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be same length for distance calculation")
        
        # Count differences and transition/transversion
        p = 0
        transitions = 0
        transversions = 0
        
        for a, b in zip(seq1, seq2):
            if a != b:
                p += 1
                if self._is_transition(a, b):
                    transitions += 1
                else:
                    transversions += 1
        
        p_frac = p / len(seq1)
        
        if self.model == "jukes_cantor":
            return self._jukes_cantor_distance(p_frac)
        elif self.model == "kimura2p":
            return self._kimura2p_distance(transitions, transversions, len(seq1))
        else:
            # Simple Hamming distance
            return p_frac
    
    @staticmethod
    def _is_transition(a: str, b: str) -> bool:
        """Check if substitution is a transition (purine<->purine, pyrimidine<->pyrimidine)."""
        transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
        return (a, b) in transitions
    
    @staticmethod
    def _jukes_cantor_distance(p: float) -> float:
        """Jukes-Cantor distance correction."""
        if p >= 0.75:
            return 10.0  # Undefined; return large value
        return -0.75 * np.log(1.0 - 4.0 * p / 3.0)
    
    @staticmethod
    def _kimura2p_distance(transitions: int, transversions: int, seq_len: int) -> float:
        """Kimura 2-parameter distance."""
        if seq_len == 0:
            return 0.0
        P = transitions / seq_len
        Q = transversions / seq_len
        
        if (1 - 2*P - Q) <= 0 or (1 - 2*Q) <= 0:
            return 10.0
        
        return (
            -0.5 * np.log(1 - 2*P - Q) - 0.25 * np.log(1 - 2*Q)
        )
    
    def neighbor_joining(self, distances: np.ndarray, labels: List[str]) -> str:
        """
        Construct neighbor-joining tree.
        
        Algorithm:
        1. Compute net divergence for each sequence pair
        2. Find closest pair (minimum corrected distance)
        3. Create internal node, update distances
        4. Repeat until 2 sequences remain
        
        Returns: Newick format tree string
        """
        n = len(labels)
        if n < 2:
            return labels[0] if n == 1 else ""
        
        active = list(range(n))
        tree_nodes = {i: labels[i] for i in range(n)}
        
        node_counter = n
        D = distances.copy()
        
        while len(active) > 2:
            # Compute net divergences
            r = np.zeros(len(active))
            for i, idx_i in enumerate(active):
                r[i] = sum(D[idx_i, idx_j] for idx_j in active if idx_i != idx_j)
            
            # Find closest pair
            min_dist = float('inf')
            best_i, best_j = -1, -1
            
            for i in range(len(active)):
                for j in range(i + 1, len(active)):
                    idx_i, idx_j = active[i], active[j]
                    corrected_dist = D[idx_i, idx_j] - (r[i] + r[j]) / (len(active) - 2)
                    
                    if corrected_dist < min_dist:
                        min_dist = corrected_dist
                        best_i, best_j = i, j
            
            # Create new node
            idx_i, idx_j = active[best_i], active[best_j]
            new_label = f"({tree_nodes[idx_i]}:{D[idx_i, idx_j]/2:.6f},{tree_nodes[idx_j]}:{D[idx_i, idx_j]/2:.6f})"
            tree_nodes[node_counter] = new_label
            
            # Update distance matrix
            new_row = np.zeros(len(D))
            for k, idx_k in enumerate(active):
                if idx_k != idx_i and idx_k != idx_j:
                    new_row[idx_k] = (D[idx_i, idx_k] + D[idx_j, idx_k] - D[idx_i, idx_j]) / 2
            
            D[node_counter, :] = new_row
            D[:, node_counter] = new_row
            
            active.remove(idx_j)
            active[active.index(idx_i)] = node_counter
            node_counter += 1
        
        # Final pair
        if len(active) == 2:
            idx_i, idx_j = active[0], active[1]
            tree = f"({tree_nodes[idx_i]}:{D[idx_i, idx_j]/2:.6f},{tree_nodes[idx_j]}:{D[idx_i, idx_j]/2:.6f});"
        else:
            tree = tree_nodes[active[0]] + ";"
        
        return tree
    
    def write_newick(self, tree_string: str, output_file: str):
        """Write tree in Newick format."""
        with open(output_file, 'w') as f:
            f.write(tree_string)
