import gzip
import mmap
from collections import defaultdict
from vcf_consensus.logger import logger

class VCFParser:
    """Parses a VCF file and retrieves variants efficiently."""

    def __init__(self, vcf_path, fasta_chromosomes, chrom_map=None):
        """
        Initializes the VCFParser and creates an index for fast access.

        Args:
            vcf_path (str): Path to the VCF file.
            fasta_chromosomes (set): Chromosome names from the FASTA reference.
            chrom_map (dict, optional): Mapping of VCF chromosome names to FASTA names.
        """
        self.vcf_path = vcf_path
        self.fasta_chromosomes = fasta_chromosomes
        self.chrom_map = chrom_map or {}
        self.sample_names = []
        self.index = defaultdict(list)  

        self._build_index()
        self.mmap_file = None
        self._open_mmap()

    def _open_file(self):
        """Opens a VCF file."""
        return gzip.open(self.vcf_path, "rt", encoding="utf-8") if self.vcf_path.endswith(".gz") else open(self.vcf_path, "r", encoding="utf-8")

    def _open_mmap(self):
        """Opens VCF with memory-mapped file for faster access."""
        with open(self.vcf_path, "r+b") as f:
            self.mmap_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

    def _build_index(self):
        """Indexes VCF positions."""
        logger.info(f"Indexing VCF: {self.vcf_path}")

        with self._open_file() as f:
            file_offset = 0  
            for line in f:
                if line.startswith("##"):
                    file_offset += len(line)
                    continue  

                if line.startswith("#CHROM"):
                    self.sample_names = line.strip().split("\t")[9:]
                    file_offset += len(line)
                    continue  

                fields = line.strip().split("\t")
                chrom, pos = fields[0], int(fields[1])

                if chrom in self.chrom_map:
                    chrom = self.chrom_map[chrom]

                if chrom not in self.fasta_chromosomes:
                    file_offset += len(line)
                    continue  

                self.index[chrom].append((pos, file_offset))
                file_offset += len(line)

        logger.info(f"Indexed {sum(len(v) for v in self.index.values())} variants from VCF.")

    def get_variants(self, chrom, start=None, end=None):
        """Retrieves variants for a given chromosome."""
        if chrom not in self.index:
            return {}

        result = {}
        for pos, offset in self.index[chrom]:
            if start is not None and pos < start:
                continue
            if end is not None and pos > end:
                break  

            variant = self._read_variant(offset)
            if variant:
                result[pos] = variant

        return result

    def _read_variant(self, offset):
        """Reads a single variant record from VCF."""
        self.mmap_file.seek(offset)
        line = self.mmap_file.readline().decode().strip()
        
        fields = line.split("\t")
        if len(fields) < 10:
            return None  

        chrom, pos, _, ref, alt, _, _, _, format_info, *samples = fields
        pos = int(pos)
        alts = alt.split(",")
        
        return {
            "REF": ref,
            "ALT": alts,
            "samples": {self.sample_names[i]: gt for i, gt in enumerate(samples)}
        }

    def get_sample_names(self):
        """Returns the list of sample names from the VCF file."""
        return self.sample_names