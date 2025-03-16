import random
from vcf_consensus.fasta_parser import FASTAParser
from vcf_consensus.vcf_parser import VCFParser
from vcf_consensus.logger import logger


def load_data(vcf_path, fasta_path, chrom_map):
    """Loads FASTA and VCF data using object-oriented parsing.

    Args:
        vcf_path (str): Path to the VCF file.
        fasta_path (str): Path to the FASTA reference genome.
        chrom_map (dict, optional): Mapping of VCF chromosome names to FASTA names.

    Returns:
        tuple: (FASTAParser instance, VCFParser instance)
    """
    fasta_parser = FASTAParser(fasta_path)
    vcf_parser = VCFParser(vcf_path, fasta_parser.get_chromosomes(), chrom_map)
    return fasta_parser, vcf_parser


def generate_single_consensus(fasta_parser, vcf_parser, length, threshold):
    """Generates a single consensus sequence.

    Args:
        fasta_parser (FASTAParser): Instance of FASTAParser.
        vcf_parser (VCFParser): Instance of VCFParser.
        length (int): Length of the consensus sequence.
        threshold (float): Allele frequency threshold.

    Returns:
        str: Formatted consensus sequence in FASTA format.
    """
    chrom = random.choice(list(fasta_parser.get_chromosomes()))
    start = random.randint(0, len(fasta_parser.get_sequence(chrom, 0, None)) - length)

    consensus = list(fasta_parser.get_sequence(chrom, start, length))

    # Выбираем случайный образец для всей последовательности
    vcf_data = vcf_parser.get_variants()
    if chrom not in vcf_data:
        return f">consensus_{chrom}_{start}\n{''.join(consensus)}"

    all_samples = list(vcf_data[chrom].values())[0]["samples"]
    sample_id = random.choice(range(len(all_samples)))  

    # Кэшируем позиции VCF
    variants_in_range = {pos: vcf_data[chrom][pos] for pos in range(start, start + length) if pos in vcf_data[chrom]}

    for pos, variant in variants_in_range.items():
        alt_alleles = variant["ALT"]
        ref_allele = variant["REF"]

        # Проверяем, несёт ли выбранный образец вариант
        genotype = variant["samples"][sample_id]
        alleles = genotype.split("/")

        selected_allele = ref_allele
        if "1" in alleles and random.random() < threshold:
            selected_allele = random.choice(alt_alleles)

        consensus[pos - start] = selected_allele  # Быстрое изменение списка

    return f">consensus_{chrom}_{start}\n{''.join(consensus)}"


def generate_consensus_sequences(vcf_path, fasta_path, length, count, threshold, output_path, seed=None, chrom_map=None):
    """Generates multiple consensus sequences from VCF and FASTA.

    Args:
        vcf_path (str): Path to the VCF file.
        fasta_path (str): Path to the FASTA reference genome.
        length (int): Length of each consensus sequence.
        count (int): Number of consensus sequences to generate.
        threshold (float): Allele frequency threshold for variant inclusion.
        output_path (str): Path to the output FASTA file.
        seed (int, optional): Random seed for reproducibility.
        chrom_map (dict, optional): Mapping of VCF chromosome names to FASTA chromosome names.
    """
    if seed:
        random.seed(seed)

    logger.info("Starting consensus sequence generation")

    fasta_parser, vcf_parser = load_data(vcf_path, fasta_path, chrom_map)
    output_sequences = []

    for i in range(count):
        consensus = generate_single_consensus(fasta_parser, vcf_parser, length, threshold)
        output_sequences.append(consensus)
        logger.info(f"Generated consensus {i + 1}/{count}")

    # Записываем результат в файл
    with open(output_path, "w") as f:
        f.write("\n".join(output_sequences) + "\n")

    logger.info(f"Saved {count} consensus sequences to {output_path}")