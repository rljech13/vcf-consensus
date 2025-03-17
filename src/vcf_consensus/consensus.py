import random
from collections import defaultdict
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


def get_consensus_positions(fasta_parser, chrom, length, count, mode="random"):
    """Generates start positions for consensus sequences.

    Ensures unique positions and limits count to available positions in `sequential` mode.

    Args:
        fasta_parser (FASTAParser): Object for FASTA access.
        chrom (str): Chromosome name.
        length (int): Length of consensus sequences.
        count (int): Number of consensus sequences.
        mode (str): "random" (random starts) or "sequential" (ordered starts).

    Returns:
        list: List of unique start positions.
    """
    chrom_length = len(fasta_parser.get_sequence(chrom, 0, None))

    if chrom_length < length:
        logger.warning(f"Chromosome {chrom} is too short for the requested consensus length ({chrom_length} < {length}). Skipping.")
        return []

    positions = []
    step = length // 2  # Default to 50% overlap

    if mode == "random":
        attempts = 0
        while len(positions) < count:
            start = random.randint(0, chrom_length - length)
            if start not in positions:
                positions.append(start)
            attempts += 1
            if attempts > count * 10:
                logger.warning(f"Could not generate enough unique positions for {chrom}. Using {len(positions)} instead of {count}.")
                break

    elif mode == "sequential":
        positions = list(range(0, chrom_length - length + 1, step))
        if len(positions) < count:
            logger.warning(f"Only {len(positions)} positions available for {chrom}, requested {count}. Using {len(positions)}.")

        positions = positions[:min(count, len(positions))]

    return sorted(positions)


def generate_single_consensus(fasta_parser, vcf_parser, chrom, start, length, threshold):
    """Generates a single consensus sequence, handling SNPs and indels.

    Args:
        fasta_parser (FASTAParser): FASTA file parser.
        vcf_parser (VCFParser): VCF file parser.
        chrom (str): Chromosome name.
        start (int): Start position.
        length (int): Length of the consensus sequence.
        threshold (float): Allele frequency threshold.

    Returns:
        str: Formatted consensus sequence in FASTA format.
    """
    consensus = list(fasta_parser.get_sequence(chrom, start, length))

    vcf_data = vcf_parser.get_variants()
    if chrom not in vcf_data:
        return f">consensus_{chrom}_{start}\n{''.join(consensus)}"

    sample_variants = [variant["samples"] for variant in vcf_data[chrom].values()]
    sample_names = list(set().union(*sample_variants))

    if not sample_names:
        logger.warning(f"No samples found for {chrom}. Using reference sequence.")
        return f">consensus_{chrom}_{start}\n{''.join(consensus)}"

    selected_sample = random.choice(sample_names)
    variants_in_range = {pos: vcf_data[chrom][pos] for pos in range(start, start + length) if pos in vcf_data[chrom]}

    shift = 0  # Track insertions and deletions

    for pos, variant in sorted(variants_in_range.items()):
        ref_allele = variant["REF"]
        alt_alleles = variant["ALT"]
        sample_data = variant["samples"]

        if selected_sample not in sample_data:
            continue

        genotype = sample_data[selected_sample]
        alleles = genotype.split("/")

        selected_allele = ref_allele
        if "1" in alleles and random.random() < threshold:
            selected_allele = random.choice(alt_alleles)

        index = pos - start + shift
        if not (0 <= index < len(consensus)):
            continue

        if len(selected_allele) > len(ref_allele):  # Insertion
            consensus[index] = selected_allele
            shift += len(selected_allele) - len(ref_allele)

        elif len(selected_allele) < len(ref_allele):  # Deletion
            del consensus[index : index + len(ref_allele)]
            shift -= len(ref_allele) - len(selected_allele)

        else:  # SNP
            consensus[index] = selected_allele

    return f">consensus_{chrom}_{start}\n{''.join(consensus)}"


def generate_consensus_sequences(
    vcf_path, fasta_path, length, count, threshold, output_path, seed=None, chrom_map=None, mode="random"
):
    """Generates consensus sequences from VCF and FASTA.

    Args:
        vcf_path (str): Path to the VCF file.
        fasta_path (str): Path to the FASTA reference genome.
        length (int): Length of each consensus sequence.
        count (int): Number of consensus sequences to generate.
        threshold (float): Allele frequency threshold for variant inclusion.
        output_path (str): Path to the output FASTA file.
        seed (int, optional): Random seed for reproducibility.
        chrom_map (dict, optional): Mapping of VCF chromosome names to FASTA chromosome names.
        mode (str, optional): "random" or "sequential" (consensus start mode).
    """
    if seed:
        random.seed(seed)

    logger.info(f"Starting consensus sequence generation in {mode} mode.")

    fasta_parser, vcf_parser = load_data(vcf_path, fasta_path, chrom_map)
    output_sequences = []

    valid_chromosomes = list(fasta_parser.get_chromosomes() & vcf_parser.get_variants().keys())

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