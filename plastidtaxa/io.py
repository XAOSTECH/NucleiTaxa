"""
Sequence I/O Module

Reading and writing sequence files (FASTQ, FASTA) with metadata tracking.
Handles both gzipped and uncompressed formats.
"""

import gzip
from pathlib import Path
from typing import Iterator, Tuple, Optional
from Bio.SeqIO.QualityIO import FastqGeneralIterator


class FastqRecord:
    """Simple FASTQ record container."""
    def __init__(self, id: str, seq: str, qual: str):
        self.id = id
        self.seq = seq
        self.qual = qual
        self.length = len(seq)
    
    def __repr__(self):
        return f"FastqRecord({self.id}, len={self.length})"


class FastqReader:
    """
    Efficient FASTQ reader with quality metrics.
    
    Supports gzip-compressed files. Yields records one at a time to minimize memory.
    
    Example:
        reader = FastqReader("reads.fastq.gz")
        for record in reader:
            print(record.id, len(record.seq), record.qual[:10])
    """
    
    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        self.is_gzipped = self.filepath.suffix == ".gz"
    
    def __iter__(self) -> Iterator[FastqRecord]:
        """Iterate over FASTQ records."""
        if self.is_gzipped:
            handle = gzip.open(self.filepath, 'rt')
        else:
            handle = open(self.filepath, 'r')
        
        try:
            for title, seq, qual in FastqGeneralIterator(handle):
                # Extract sequence ID (remove trailing /1 or /2 from paired-end reads)
                seq_id = title.split()[0] if title else title
                yield FastqRecord(seq_id, seq, qual)
        finally:
            handle.close()
    
    def count_reads(self) -> int:
        """Count total reads in file (expensive for large files)."""
        count = 0
        for _ in self:
            count += 1
        return count
    
    def quality_stats(self) -> dict:
        """
        Compute quality score statistics across all reads.
        
        Returns: dict with keys: mean_phred, min_phred, max_phred, n_reads
        """
        phred_scores = []
        n_reads = 0
        
        for record in self:
            n_reads += 1
            # Convert ASCII quality to Phred scores
            scores = [ord(c) - 33 for c in record.qual]
            phred_scores.extend(scores)
        
        if not phred_scores:
            return {"mean_phred": 0, "min_phred": 0, "max_phred": 0, "n_reads": 0}
        
        import statistics
        return {
            "mean_phred": statistics.mean(phred_scores),
            "min_phred": min(phred_scores),
            "max_phred": max(phred_scores),
            "n_reads": n_reads,
        }


class FastaWriter:
    """
    FASTA writer for sequences.
    
    Example:
        writer = FastaWriter("output.fasta")
        writer.write("seq1", "ACGTACGT")
        writer.write("seq2", "TGCATGCA")
        writer.close()
    """
    
    def __init__(self, filepath: str, gzip_output: bool = False):
        self.filepath = Path(filepath)
        self.gzip_output = gzip_output
        
        if gzip_output:
            self.handle = gzip.open(self.filepath, 'wt')
        else:
            self.handle = open(self.filepath, 'w')
    
    def write(self, seq_id: str, sequence: str, description: str = ""):
        """Write a FASTA record."""
        header = f">{seq_id}"
        if description:
            header += f" {description}"
        self.handle.write(header + "\n")
        
        # Wrap sequence at 80 characters
        for i in range(0, len(sequence), 80):
            self.handle.write(sequence[i:i+80] + "\n")
    
    def close(self):
        """Close the output file."""
        self.handle.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        self.close()


class FastqFilter:
    """
    FASTQ filtering with configurable quality control.
    
    Bioinformatic rationale:
    - Low-quality reads introduce errors in taxonomic inference
    - Ambiguous bases (N) complicate alignment and classification
    - Trimming removes adapters and degraded sequence ends
    - Expected error filtering is more sensitive than simple Q-score cutoffs
    """
    
    def __init__(self, trim_left: Tuple[int, int] = (0, 0),
                 trim_right: Tuple[int, int] = (0, 0),
                 trunc_len: Tuple[int, int] = (0, 0),
                 max_ee: Tuple[float, float] = (2.0, 2.0),
                 min_quality: int = 3,
                 max_n: int = 0):
        """
        Args:
            trim_left: (fwd, rev) nucleotides to trim from 5' end
            trim_right: (fwd, rev) nucleotides to trim from 3' end
            trunc_len: (fwd, rev) hard truncation position (0=disable)
            max_ee: (fwd, rev) max expected error threshold
            min_quality: minimum mean quality score
            max_n: max ambiguous bases allowed
        """
        self.trim_left = trim_left
        self.trim_right = trim_right
        self.trunc_len = trunc_len
        self.max_ee = max_ee
        self.min_quality = min_quality
        self.max_n = max_n
    
    def calculate_expected_error(self, qual: str) -> float:
        """
        Calculate expected error (EE) for a sequence.
        
        Bioinformatic principle: EE is the sum of error probabilities.
        Lower EE = higher confidence in sequence accuracy.
        
        Formula: EE = sum(10^(-Q/10)) for each base
        """
        ee = 0.0
        for q_char in qual:
            q = ord(q_char) - 33  # Convert ASCII to Phred
            ee += 10 ** (-q / 10)
        return ee
    
    def filter_record(self, record: FastqRecord, is_forward: bool = True) -> Optional[FastqRecord]:
        """
        Filter a single FASTQ record.
        
        Returns None if record fails QC, else returns filtered record.
        """
        seq = record.seq
        qual = record.qual
        
        # Apply trimming
        trim_l, trim_r = self.trim_left if is_forward else self.trim_left
        if trim_l or trim_r:
            seq = seq[trim_l:len(seq)-trim_r if trim_r else None]
            qual = qual[trim_l:len(qual)-trim_r if trim_r else None]
        
        # Apply truncation
        trunc = self.trunc_len[0] if is_forward else self.trunc_len[1]
        if trunc > 0 and len(seq) > trunc:
            seq = seq[:trunc]
            qual = qual[:trunc]
        
        # Check minimum length
        if len(seq) < 10:  # Arbitrary minimum
            return None
        
        # Check max N's
        if seq.count('N') > self.max_n:
            return None
        
        # Check expected error
        max_ee = self.max_ee[0] if is_forward else self.max_ee[1]
        ee = self.calculate_expected_error(qual)
        if ee > max_ee:
            return None
        
        # Check mean quality
        mean_q = sum(ord(c) - 33 for c in qual) / len(qual)
        if mean_q < self.min_quality:
            return None
        
        return FastqRecord(record.id, seq, qual)
