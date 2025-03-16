import gzip
from vcf_consensus.logger import logger

class FASTAParser:
    """Parses a FASTA file (supports .fasta and .fasta.gz) and stores sequences in memory."""

    def __init__(self, fasta_path):
        """
        Initializes the FASTAParser and loads sequences.

        Args:
            fasta_path (str): Path to the FASTA file (supports .fasta and .fasta.gz).
        """
        self.fasta_path = fasta_path
        self.sequences = {}
        self.fasta_chromosomes = set()
        self._parse_fasta()

    def _open_file(self):
        """Opens a FASTA file, handling both .fasta and .fasta.gz formats."""
        if self.fasta_path.endswith(".gz"):
            return gzip.open(self.fasta_path, "rt", encoding="utf-8")
        return open(self.fasta_path, "r", encoding="utf-8")

    def _parse_fasta(self):
        """Internal method to parse the FASTA file and store sequences."""
        logger.info(f"Loading FASTA: {self.fasta_path}")
        current_chrom = None
        current_seq = []

        with self._open_file() as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_chrom:
                        self.sequences[current_chrom] = "".join(current_seq)
                    current_chrom = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_chrom:
            self.sequences[current_chrom] = "".join(current_seq)

        self.fasta_chromosomes = set(self.sequences.keys())
        logger.info(f"Loaded {len(self.fasta_chromosomes)} chromosomes from FASTA")

    def get_sequence(self, chrom, start=0, length=None):
        """
        Retrieves a sequence from the FASTA file.

        Args:
            chrom (str): Chromosome name.
            start (int): Start position.
            length (int, optional): Length of sequence. If None, returns the full chromosome.

        Returns:
            str: The extracted sequence.
        """
        sequence = self.sequences.get(chrom, "")
        if length is None:
            return sequence[start:]
        return sequence[start:start + length]

    def get_chromosomes(self):
        """Returns a set of chromosome names from the FASTA file."""
        return self.fasta_chromosomes