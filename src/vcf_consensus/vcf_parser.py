import gzip
from collections import defaultdict
from vcf_consensus.logger import logger

class VCFParser:
    """Parses a VCF file and efficiently retrieves variant data."""

    def __init__(self, vcf_path, fasta_chromosomes, chrom_map=None):
        """
        Initializes the VCFParser and loads relevant variant data.

        Args:
            vcf_path (str): Path to the VCF file.
            fasta_chromosomes (set): Chromosome names from the FASTA reference.
            chrom_map (dict, optional): Mapping of VCF chromosome names to FASTA names.
        """
        self.vcf_path = vcf_path
        self.fasta_chromosomes = fasta_chromosomes
        self.chrom_map = chrom_map or {}
        self.vcf_data = defaultdict(dict)
        self.sample_names = []

        self._parse_vcf()

    def _open_file(self):
        """Opens a VCF file, handling both .vcf and .vcf.gz formats."""
        return gzip.open(self.vcf_path, "rt", encoding="utf-8") if self.vcf_path.endswith(".gz") else open(self.vcf_path, "r", encoding="utf-8")

    def _parse_vcf(self):
        """Parses the VCF file and stores only necessary variant data."""
        logger.info(f"Loading VCF: {self.vcf_path}")

        with self._open_file() as f:
            for line in f:
                if line.startswith("##"):
                    continue

                if line.startswith("#CHROM"):
                    self.sample_names = line.strip().split("\t")[9:]
                    continue

                fields = line.strip().split("\t")
                chrom, pos, _, ref, alt, _, _, _, format_info, *samples = fields
                pos = int(pos)
                alts = alt.split(",")

                if chrom in self.chrom_map:
                    chrom = self.chrom_map[chrom]

                if chrom not in self.fasta_chromosomes:
                    continue  

                self.vcf_data[chrom][pos] = {
                    "REF": ref,
                    "ALT": alts,
                    "samples": {self.sample_names[i]: gt for i, gt in enumerate(samples)}
                }

        logger.info(f"Loaded {sum(len(v) for v in self.vcf_data.values())} variants from VCF.")

    def get_variants(self):
        """Returns the parsed VCF data."""
        return self.vcf_data

    def get_sample_names(self):
        """Returns the list of sample names from the VCF file."""
        return self.sample_names