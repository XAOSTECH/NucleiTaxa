"""
Safe Sequence Fetching Module

Retrieves sequences from remote databases (NCBI Entrez, SILVA, etc.)
with XXE (XML External Entity) injection protection.

Security Note: Uses defusedxml to parse XML safely, mitigating
the XXE vulnerability in Biopython <= 1.86.
"""

from typing import List, Dict, Optional
from defusedxml import ElementTree as ET
from io import StringIO
import time

try:
    from Bio import Entrez
except ImportError:
    Entrez = None


def safe_entrez_fetch(query: str, db: str = "nucleotide",
                      retmax: int = 100, delay: float = 0.5) -> List[Dict]:
    """
    Safely fetch sequences from NCBI Entrez using XXE-protected parsing.
    
    Bioinformatic use case:
    - Retrieve reference sequences for taxonomic comparison
    - Download SILVA/RDP training datasets
    - Query reference genomes for validation
    
    Security: Uses defusedxml to prevent XXE injection attacks that
    were possible in Biopython <= 1.86 via Bio.Entrez.read()
    
    Args:
        query: Entrez search query (e.g., "16S ribosomal RNA[Title]")
        db: Database ("nucleotide", "protein", "taxonomy")
        retmax: Max results to return
        delay: Delay between requests (seconds) to respect NCBI rate limits
    
    Returns: List of dicts with 'id', 'title', 'sequence' keys
    """
    if Entrez is None:
        raise ImportError("Biopython not installed. Install with: pip install biopython")
    
    # Set email for NCBI (required by their guidelines)
    Entrez.email = "bioinformatics@example.com"
    
    results = []
    
    try:
        # Search phase: Get UIDs matching query
        print(f"Searching NCBI {db} for: {query}")
        search_handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        
        # Use defusedxml to parse search results safely
        search_xml = search_handle.read()
        search_root = ET.fromstring(search_xml)
        
        uids = []
        for uid_elem in search_root.findall('.//Id'):
            if uid_elem.text:
                uids.append(uid_elem.text)
        
        print(f"Found {len(uids)} sequences")
        
        # Fetch phase: Download full records
        if not uids:
            return results
        
        time.sleep(delay)  # Rate limiting
        fetch_handle = Entrez.efetch(db=db, id=",".join(uids), rettype="gb", retmode="xml")
        fetch_xml = fetch_handle.read()
        
        # Parse with defusedxml
        fetch_root = ET.fromstring(fetch_xml)
        
        # Extract sequence records (structure varies by db)
        for seq_elem in fetch_root.findall('.//GBSeq') or fetch_root.findall('.//Seq'):
            record = _parse_genbank_record(seq_elem)
            if record:
                results.append(record)
        
        return results
    
    except Exception as e:
        print(f"Error fetching sequences: {e}")
        return []


def _parse_genbank_record(elem) -> Optional[Dict]:
    """Parse GenBank XML element safely."""
    try:
        record = {}
        
        # Extract metadata
        locus_elem = elem.find('.//GBSeq_locus')
        if locus_elem is not None:
            record['id'] = locus_elem.text
        
        definition_elem = elem.find('.//GBSeq_definition')
        if definition_elem is not None:
            record['title'] = definition_elem.text
        
        # Extract sequence
        seq_elem = elem.find('.//GBSeq_sequence')
        if seq_elem is not None:
            record['sequence'] = seq_elem.text.upper()
        
        return record if 'id' in record else None
    
    except Exception:
        return None


def fetch_silva_reference(target: str = "23S", version: str = "138.1") -> str:
    """
    Download SILVA reference database for taxonomic assignment.
    
    Args:
        target: rRNA target ("23S", "16S", "18S")
        version: SILVA database version
    
    Returns: Path to downloaded FASTA file
    
    Note: Large files (>1GB) - consider pre-downloading or using local copies
    """
    import urllib.request
    from pathlib import Path
    
    silva_urls = {
        "23S": f"https://www.arb-silva.de/fileadmin/silva_databases/release_{version}/Exports/SILVA_{version}_LSURef_tax_silva.fasta.gz",
        "16S": f"https://www.arb-silva.de/fileadmin/silva_databases/release_{version}/Exports/SILVA_{version}_SSURef_tax_silva.fasta.gz",
    }
    
    if target not in silva_urls:
        raise ValueError(f"Unknown target: {target}. Available: {list(silva_urls.keys())}")
    
    url = silva_urls[target]
    output_path = Path(f"silva_{target}_{version}.fasta.gz")
    
    print(f"Downloading {target} from SILVA v{version}...")
    print(f"URL: {url}")
    print("This may take several minutes...")
    
    urllib.request.urlretrieve(url, output_path)
    print(f"Downloaded to: {output_path}")
    
    return str(output_path)
