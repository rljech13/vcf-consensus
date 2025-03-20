import mmap
import gzip
from vcf_consensus.logger import logger

class FASTAParser:
    """Parses a FASTA file using memory-mapped access for efficient sequence extraction."""

    def __init__(self, fasta_path):
        """
        Initializes the FASTAParser and maps the FASTA file into memory.

        Args:
            fasta_path (str): Path to the FASTA file.
        """
        self.fasta_path = fasta_path
        self.index = {}
        self.mmap_file = None
        self._parse_fasta_headers()

    def _open_file(self):
        """Opens a FASTA file (supports .fasta and .fasta.gz)."""
        return gzip.open(self.fasta_path, "rt", encoding="utf-8") if self.fasta_path.endswith(".gz") else open(self.fasta_path, "r", encoding="utf-8")

    def _parse_fasta_headers(self):
        """Indexes FASTA headers and chromosome lengths using memory-mapped access."""
        logger.info(f"Indexing FASTA: {self.fasta_path}")

        with open(self.fasta_path, "r", encoding="utf-8") as f:
            self.mmap_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            chrom = None
            seq_start = 0
            seq_length = 0

            for line in iter(self.mmap_file.readline, b""):
                line = line.decode().strip()

                if line.startswith(">"):
                    if chrom:
                        self.index[chrom]["length"] = seq_length
                    chrom = line[1:].split()[0]
                    self.index[chrom] = {"position": seq_start, "length": 0}
                    seq_start = self.mmap_file.tell()
                    seq_length = 0
                else:
                    seq_length += len(line)

            if chrom:
                self.index[chrom]["length"] = seq_length

        logger.info(f"Indexed {len(self.index)} chromosomes from FASTA.")

    def get_sequence(self, chrom, start=0, length=None):
        """Retrieves a sequence from the memory-mapped FASTA file."""
        if chrom not in self.index:
            raise ValueError(f"Chromosome {chrom} not found in FASTA!")

        chrom_position = self.index[chrom]["position"]
        chrom_length = self.index[chrom]["length"]

        if start >= chrom_length:
            raise ValueError(f"Start position {start} is out of bounds for chromosome {chrom}")

        if length is None or start + length > chrom_length:
            length = chrom_length - start

        self.mmap_file.seek(chrom_position + start)
        sequence = self.mmap_file.read(length).decode().replace("\n", "")

        return sequence[:length]

    def get_chromosomes(self):
        """Returns indexed chromosome names."""
        return set(self.index.keys())

    def get_chromosome_length(self, chrom):
        """Returns the length of the specified chromosome."""
        if chrom not in self.index:
            raise ValueError(f"Chromosome {chrom} not found in FASTA!")
        return self.index[chrom]["length"]