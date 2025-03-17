import gzip
from vcf_consensus.logger import logger

class FASTAParser:
    """Parses a FASTA file (supports .fasta and .fasta.gz) efficiently."""

    def __init__(self, fasta_path):
        """
        Initializes the FASTAParser and indexes chromosome positions.

        Args:
            fasta_path (str): Path to the FASTA file (supports .fasta and .fasta.gz).
        """
        self.fasta_path = fasta_path
        self.index = {}  # Stores chromosome start positions
        self._parse_fasta_headers()

    def _open_file(self):
        """Opens a FASTA file, handling both .fasta and .fasta.gz formats."""
        return gzip.open(self.fasta_path, "rt", encoding="utf-8") if self.fasta_path.endswith(".gz") else open(self.fasta_path, "r", encoding="utf-8")

    def _parse_fasta_headers(self):
        """Indexes FASTA headers to enable efficient sequence extraction."""
        logger.info(f"Indexing FASTA: {self.fasta_path}")

        with self._open_file() as f:
            position = 0
            chrom = None
            for line in f:
                if line.startswith(">"):
                    chrom = line[1:].split()[0]
                    self.index[chrom] = position
                position += len(line)

        logger.info(f"Indexed {len(self.index)} chromosomes from FASTA.")

    def get_sequence(self, chrom, start=0, length=None):
        """Extracts a sequence from the FASTA file efficiently.

        Args:
            chrom (str): Chromosome name.
            start (int): Start position.
            length (int, optional): Length of sequence.

        Returns:
            str: Extracted DNA sequence.
        """
        if chrom not in self.index:
            raise ValueError(f"Chromosome {chrom} not found in FASTA!")

        sequence = []
        reading = False

        with self._open_file() as f:
            for line in f:
                if line.startswith(">"):
                    reading = line[1:].split()[0] == chrom
                    continue
                if reading:
                    sequence.append(line.strip())
                    if length and len("".join(sequence)) >= start + length:
                        break

        full_seq = "".join(sequence)
        return full_seq[start:start + length] if length else full_seq

    def get_chromosomes(self):
        """Returns the indexed chromosome names."""
        return set(self.index.keys())