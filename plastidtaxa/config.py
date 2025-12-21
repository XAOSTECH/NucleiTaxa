"""
Nucleid Profile Configuration

Defines target sequences (plastid 23S/16S, mitochondrial, ITS, etc.)
with database paths, quality thresholds, and processing parameters.
"""

from dataclasses import dataclass
from typing import Optional, Dict, List


@dataclass
class NuclidProfile:
    """
    Configuration profile for a nucleid target.
    
    Attributes:
        name: Profile identifier (e.g., "plastid_23s", "mito_12s", "fungal_its")
        description: Human-readable description
        region: rRNA region or target (e.g., "23S", "16S", "ITS2")
        expected_length: Expected sequence length in bp
        trim_left: Nucleotides to trim from left (primer/adapter removal)
        trim_right: Nucleotides to trim from right
        trunc_len_f: Truncate forward reads at this position
        trunc_len_r: Truncate reverse reads at this position
        max_ee: Max expected error threshold
        min_quality: Minimum mean quality score
        max_n: Maximum ambiguous nucleotides allowed
        reference_db: Path to reference database (FASTA)
        silva_train: Path to SILVA training set (for RDP classifier)
    """
    name: str
    description: str
    region: str
    expected_length: int
    trim_left: tuple = (0, 0)  # (forward, reverse)
    trim_right: tuple = (0, 0)
    trunc_len_f: int = 0  # 0 = no truncation
    trunc_len_r: int = 0
    max_ee: tuple = (2.0, 2.0)  # (forward, reverse)
    min_quality: int = 3
    max_n: int = 0
    reference_db: Optional[str] = None
    silva_train: Optional[str] = None
    
    def __str__(self) -> str:
        return f"{self.name} ({self.region}): {self.description}"


# Pre-defined profiles
PROFILES: Dict[str, NuclidProfile] = {
    "plastid_23s": NuclidProfile(
        name="plastid_23s",
        description="Plastid 23S ribosomal RNA (chloroplast large subunit)",
        region="23S",
        expected_length=2900,
        trim_left=(21, 20),
        trim_right=(0, 0),
        trunc_len_f=285,
        trunc_len_r=220,
        max_ee=(2.0, 2.0),
        min_quality=3,
        max_n=0,
    ),
    "plastid_16s": NuclidProfile(
        name="plastid_16s",
        description="Plastid 16S ribosomal RNA (chloroplast small subunit)",
        region="16S",
        expected_length=1500,
        trim_left=(15, 15),
        trim_right=(0, 0),
        trunc_len_f=250,
        trunc_len_r=200,
        max_ee=(2.0, 2.0),
        min_quality=3,
        max_n=0,
    ),
    "mito_12s": NuclidProfile(
        name="mito_12s",
        description="Mitochondrial 12S ribosomal RNA (animal mtDNA)",
        region="12S",
        expected_length=950,
        trim_left=(10, 10),
        trim_right=(0, 0),
        trunc_len_f=200,
        trunc_len_r=180,
        max_ee=(2.0, 2.0),
        min_quality=3,
        max_n=0,
    ),
    "fungal_its": NuclidProfile(
        name="fungal_its",
        description="Fungal Internal Transcribed Spacer (ITS1/ITS2)",
        region="ITS2",
        expected_length=300,
        trim_left=(0, 0),
        trim_right=(0, 0),
        trunc_len_f=0,  # Variable length region
        trunc_len_r=0,
        max_ee=(3.0, 3.0),
        min_quality=3,
        max_n=1,
    ),
    "bacterial_16s": NuclidProfile(
        name="bacterial_16s",
        description="Bacterial/Archaeal 16S ribosomal RNA",
        region="16S",
        expected_length=1500,
        trim_left=(13, 12),
        trim_right=(0, 0),
        trunc_len_f=240,
        trunc_len_r=160,
        max_ee=(2.0, 2.0),
        min_quality=3,
        max_n=0,
    ),
}


def get_profile(name: str) -> NuclidProfile:
    """Retrieve a pre-defined nucleid profile by name."""
    if name not in PROFILES:
        available = ", ".join(PROFILES.keys())
        raise ValueError(f"Unknown profile '{name}'. Available: {available}")
    return PROFILES[name]


def list_profiles() -> List[str]:
    """Return list of available profile names."""
    return list(PROFILES.keys())
