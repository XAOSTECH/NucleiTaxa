"""
Command-line interface for PlastidTaxa pipeline.

Example usage:
    plastidtaxa --profile plastid_23s \\
        --forward reads_R1.fastq.gz \\
        --reverse reads_R2.fastq.gz \\
        --output ./results \\
        --reference silva_23s.fasta
"""

import click
import sys
from pathlib import Path
from .config import list_profiles, get_profile
from .io import FastqReader, FastqFilter
from .qc import QualityController
from .denoise import Denoiser
from .taxonomy import TaxonomyAssigner


@click.group()
def main():
    """PlastidTaxa: Nucleid-agnostic sequence taxonomy pipeline"""
    pass


@main.command()
def list_profiles_cmd():
    """List available nucleid profiles."""
    profiles = list_profiles()
    for prof in profiles:
        p = get_profile(prof)
        click.echo(f"{prof:20s} - {p.description}")


@main.command()
@click.option('--profile', default='plastid_23s', 
              help='Nucleid profile (plastid_23s, plastid_16s, mito_12s, fungal_its, bacterial_16s)')
@click.option('--forward', required=True, help='Forward reads (R1, gzipped FASTQ)')
@click.option('--reverse', required=True, help='Reverse reads (R2, gzipped FASTQ)')
@click.option('--reference', required=True, help='Reference database (FASTA)')
@click.option('--output', default='./results', help='Output directory')
@click.option('--min-quality', type=int, default=3, help='Minimum quality score')
@click.option('--max-ee', type=float, default=2.0, help='Max expected error')
def run_pipeline(profile, forward, reverse, reference, output, min_quality, max_ee):
    """
    Run complete PlastidTaxa analysis pipeline.
    
    Steps:
    1. Quality control & filtering
    2. Sequence denoising (ASV inference)
    3. Merge paired-end reads
    4. Taxonomic assignment
    5. Generate output tables
    """
    
    # Validate inputs
    forward_path = Path(forward)
    reverse_path = Path(reverse)
    reference_path = Path(reference)
    output_dir = Path(output)
    
    if not forward_path.exists():
        click.echo(f"Error: Forward reads not found: {forward}", err=True)
        sys.exit(1)
    
    if not reverse_path.exists():
        click.echo(f"Error: Reverse reads not found: {reverse}", err=True)
        sys.exit(1)
    
    if not reference_path.exists():
        click.echo(f"Error: Reference database not found: {reference}", err=True)
        sys.exit(1)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get nucleid profile
    try:
        prof = get_profile(profile)
        click.echo(f"Profile: {prof}")
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    
    # Step 1: Quality Control
    click.echo("\n[1/5] Quality Control...")
    qc = QualityController(prof)
    fwd_filt, rev_filt = qc.filter_fastq_pair(
        str(forward_path),
        str(reverse_path),
        str(output_dir / "filtered")
    )
    click.echo(f"  Input reads: {qc.stats['input']}")
    click.echo(f"  Passed QC: {qc.stats['passed_qc']}")
    
    # Step 2: Dereplication
    click.echo("\n[2/5] Dereplication...")
    derep_fwd = qc.dereplicate(fwd_filt)
    click.echo(f"  Unique sequences (forward): {len(derep_fwd)}")
    
    # Step 3: Denoising (ASV inference)
    click.echo("\n[3/5] Denoising (ASV inference)...")
    denoiser = Denoiser(error_threshold=0.01)
    sequences, abundance = denoiser.cluster_sequences(
        list(derep_fwd.keys()),
        list(derep_fwd.values()),
        min_abundance=2
    )
    click.echo(f"  ASVs identified: {len(sequences)}")
    
    # Step 4: Taxonomic Assignment
    click.echo("\n[4/5] Taxonomic Assignment...")
    try:
        assigner = TaxonomyAssigner(str(reference_path), method="kmer")
        taxa_assignments = assigner.batch_assign(sequences)
        click.echo(f"  Classified: {len([t for t in taxa_assignments if t['confidence'] > 0.8])}")
    except Exception as e:
        click.echo(f"  Warning: Taxonomy assignment failed: {e}", err=True)
        taxa_assignments = [{"taxonomy": "unclassified", "confidence": 0.0}] * len(sequences)
    
    # Step 5: Generate Output
    click.echo("\n[5/5] Generating outputs...")
    _write_results(sequences, abundance, taxa_assignments, output_dir)
    
    click.echo(f"\nResults written to: {output_dir}")
    click.echo("✓ Pipeline complete!")


@main.command()
@click.option('--input', required=True, help='Input FASTQ file')
@click.option('--profile', default='plastid_23s', help='Nucleid profile')
def qc_report(input, profile):
    """Generate quality control report."""
    prof = get_profile(profile)
    reader = FastqReader(input)
    stats = reader.quality_stats()
    
    click.echo(f"\nQuality Report: {Path(input).name}")
    click.echo(f"Profile: {prof}")
    click.echo(f"  Total reads: {stats['n_reads']}")
    click.echo(f"  Mean Phred: {stats['mean_phred']:.2f}")
    click.echo(f"  Min Phred: {stats['min_phred']}")
    click.echo(f"  Max Phred: {stats['max_phred']}")


def _write_results(sequences, abundance, taxa, output_dir: Path):
    """Write ASV table, taxa table, and FASTA."""
    
    # ASV table (OTU-like abundance matrix)
    asv_file = output_dir / "ASV_table.txt"
    with open(asv_file, 'w') as f:
        f.write("ASV_ID\tSequence\tAbundance\tTaxonomy\tConfidence\n")
        for i, (seq, abund) in enumerate(zip(sequences, abundance)):
            tax = taxa[i] if i < len(taxa) else {}
            taxonomy = tax.get('taxonomy', 'unclassified')
            confidence = tax.get('confidence', 0.0)
            f.write(f"ASV_{i+1}\t{seq}\t{abund}\t{taxonomy}\t{confidence:.3f}\n")
    
    click.echo(f"  ASV table: {asv_file}")
    
    # Representative sequences
    fasta_file = output_dir / "ASV_seqs.fasta"
    with open(fasta_file, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f">ASV_{i+1}\n{seq}\n")
    
    click.echo(f"  FASTA: {fasta_file}")


if __name__ == '__main__':
    main()
