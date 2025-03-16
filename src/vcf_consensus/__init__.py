"""Initialize vcf_consensus package."""
__version__ = "0.1.0"

from .fasta_parser import FASTAParser  
from .vcf_parser import VCFParser 
from .consensus import generate_consensus_sequences