import random
from vcf_consensus.fasta_parser import FASTAParser
from vcf_consensus.vcf_parser import VCFParser
from vcf_consensus.logger import logger


def load_data(vcf_path, fasta_path, chrom_map):
    """Loads FASTA and VCF data using optimized indexing.

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


def get_consensus_positions(fasta_parser, chrom, length, count, mode="random"):
    """Generates start positions for consensus sequences efficiently.

    Args:
        fasta_parser (FASTAParser): FASTA parser with indexed sequences.
        chrom (str): Chromosome name.
        length (int): Length of consensus sequences.
        count (int): Number of sequences to generate.
        mode (str): "random" (randomized start positions) or "sequential" (ordered).

    Returns:
        list: List of start positions.
    """
    chrom_length = fasta_parser.get_chromosome_length(chrom)

    if chrom_length < length:
        logger.warning(f"Chromosome {chrom} is too short for requested consensus length. Skipping.")
        return []

    positions = set()
    step = length // 2  # 50% overlap

    if mode == "random":
        while len(positions) < count:
            start = random.randint(0, chrom_length - length)
            positions.add(start)
            if len(positions) >= count:
                break

    elif mode == "sequential":
        positions = list(range(0, chrom_length - length + 1, step))
        positions = positions[:min(count, len(positions))]

    return sorted(positions)


def generate_single_consensus(fasta_parser, vcf_parser, chrom, start, length, threshold):
    """Generates a single consensus sequence with correct SNPs and indels.

    Args:
        fasta_parser (FASTAParser): Indexed FASTA parser.
        vcf_parser (VCFParser): Indexed VCF parser.
        chrom (str): Chromosome name.
        start (int): Start position.
        length (int): Length of the sequence.
        threshold (float): Allele frequency threshold.

    Returns:
        str: FASTA-formatted consensus sequence.
    """
    ref_sequence = bytearray(fasta_parser.get_sequence(chrom, start, length), "utf-8")

    # Используем оптимизированный `get_variants`, который сразу фильтрует по start, end
    variants_in_range = vcf_parser.get_variants(chrom, start, start + length)

    if not variants_in_range:
        return f">consensus_{chrom}_{start}\n{ref_sequence.decode('utf-8')}"

    sample_names = vcf_parser.get_sample_names()
    if not sample_names:
        logger.warning(f"No samples found in VCF for {chrom}. Keeping reference sequence.")
        return f">consensus_{chrom}_{start}\n{ref_sequence.decode('utf-8')}"

    selected_sample = random.choice(sample_names)
    shift = 0  # Корректировка позиций из-за инделов

    for pos, variant in sorted(variants_in_range.items()):
        ref_allele = variant["REF"]
        alt_alleles = variant["ALT"]
        sample_data = variant["samples"]

        if selected_sample not in sample_data:
            continue

        genotype = sample_data[selected_sample]
        alleles = genotype.split("/")

        if all(a == "0" for a in alleles):  # Гомозиготный референс (0/0)
            continue

        alt_indices = [int(a) for a in alleles if a.isdigit()]
        if not alt_indices:
            continue

        selected_allele = ref_allele
        if any(a > 0 and random.random() < threshold for a in alt_indices):
            selected_allele = alt_alleles[max(alt_indices) - 1]  # Берем наиболее частый ALT

        index = pos - start + shift
        if not (0 <= index < len(ref_sequence)):
            continue

        # Обработка SNP, инсерций и делеций
        if len(selected_allele) > len(ref_allele):  # Инсерция
            ref_sequence[index:index+len(ref_allele)] = bytearray(selected_allele, "utf-8")
            shift += len(selected_allele) - len(ref_allele)

        elif len(selected_allele) < len(ref_allele):  # Делеция
            del ref_sequence[index : index + len(ref_allele)]
            shift -= len(ref_allele) - len(selected_allele)

        else:  # SNP
            ref_sequence[index:index + len(ref_allele)] = bytearray(selected_allele, "utf-8")

    return f">consensus_{chrom}_{start}\n{ref_sequence.decode('utf-8')}"

def generate_consensus_sequences(
    vcf_path, fasta_path, length, count, threshold, output_path, seed=None, chrom_map=None, mode="random"
):
    """Generates consensus sequences from VCF and FASTA with optimized performance.

    Args:
        vcf_path (str): Path to the VCF file.
        fasta_path (str): Path to the FASTA reference genome.
        length (int): Length of consensus sequences.
        count (int): Number of sequences to generate.
        threshold (float): Allele frequency threshold.
        output_path (str): Path to the output FASTA file.
        seed (int, optional): Random seed for reproducibility.
        chrom_map (dict, optional): Chromosome name mapping.
        mode (str, optional): "random" or "sequential".
    """
    if seed:
        random.seed(seed)

    logger.info(f"Starting consensus sequence generation in {mode} mode.")

    fasta_parser, vcf_parser = load_data(vcf_path, fasta_path, chrom_map)
    output_sequences = []

    # Используем vcf_parser.index.keys(), т.к. get_variants() требует chrom, start, end
    valid_chromosomes = list(fasta_parser.get_chromosomes() & vcf_parser.index.keys())

    if not valid_chromosomes:
        logger.error("No matching chromosomes between FASTA and VCF!")
        return

    positions_dict = {
        chrom: get_consensus_positions(fasta_parser, chrom, length, count // len(valid_chromosomes), mode)
        for chrom in valid_chromosomes
    }

    consensus_count = 0
    for chrom, positions in positions_dict.items():
        for start in positions:
            if consensus_count >= count:
                break

            consensus = generate_single_consensus(fasta_parser, vcf_parser, chrom, start, length, threshold)
            if consensus:
                output_sequences.append(consensus)
                consensus_count += 1

            logger.info(f"Generated {consensus_count}/{count} consensus at {chrom}:{start}")

    with open(output_path, "w") as f:
        f.write("\n".join(output_sequences) + "\n")

    logger.info(f"Saved {len(output_sequences)} consensus sequences to {output_path}")