import gzip
from collections import defaultdict
from vcf_consensus.logger import logger

class VCFParser:
    """Parses a VCF file (supports .vcf and .vcf.gz) and stores variants in a dictionary.
    
    Attributes:
        vcf_data (dict): Parsed variant data {chrom: {position: variant_info}}.
        fasta_chromosomes (set): Set of chromosome names from the FASTA reference.
        chrom_map (dict): Optional mapping for chromosome names.
    """

    def __init__(self, vcf_path, fasta_chromosomes, chrom_map=None):
        """
        Initializes the VCFParser and loads variant data.

        Args:
            vcf_path (str): Path to the VCF file (supports .vcf and .vcf.gz).
            fasta_chromosomes (set): Chromosome names from the FASTA reference.
            chrom_map (dict, optional): Mapping of VCF chromosome names to FASTA names.
        """
        self.vcf_path = vcf_path
        self.fasta_chromosomes = fasta_chromosomes
        self.chrom_map = chrom_map or {}
        self.vcf_data = defaultdict(dict)
        self.vcf_chromosomes = set()

        self._parse_vcf()

    def _open_file(self):
        """Opens a VCF file, handling both .vcf and .vcf.gz formats."""
        if self.vcf_path.endswith(".gz"):
            return gzip.open(self.vcf_path, "rt", encoding="utf-8")
        return open(self.vcf_path, "r", encoding="utf-8")

    def _parse_vcf(self):
        """Internal method to parse the VCF file and store data in a dictionary."""
        logger.info(f"Loading VCF: {self.vcf_path}")
        with self._open_file() as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                chrom, pos, _, ref, alt, _, _, _, _, *samples = fields
                pos = int(pos)
                alts = alt.split(",")

                self.vcf_chromosomes.add(chrom)

                
                if chrom in self.chrom_map:
                    chrom = self.chrom_map[chrom]

                self.vcf_data[chrom][pos] = {
                    "REF": ref,
                    "ALT": alts,
                    "samples": samples
                }

        self._validate_chromosomes()

    def _validate_chromosomes(self):
        """Checks for mismatches between VCF and FASTA chromosome names."""
        unmatched_chroms = self.vcf_chromosomes - self.fasta_chromosomes
        if unmatched_chroms:
            logger.warning("Chromosome names in VCF do not match FASTA reference!")
            logger.warning(f"VCF chromosomes: {self.vcf_chromosomes}")
            logger.warning(f"FASTA chromosomes: {self.fasta_chromosomes}")
            logger.warning(f"Unmatched chromosomes: {unmatched_chroms}")
            logger.warning("Use --chrom-map to specify a mapping (e.g., '1=chr1,2=chr2').")

        logger.info(f"Loaded {sum(len(v) for v in self.vcf_data.values())} variants from VCF")

    def get_variants(self):
        """Returns the parsed VCF data."""
        return self.vcf_data