import argparse
import time
from vcf_consensus.consensus import generate_consensus_sequences
from vcf_consensus.logger import logger

def parse_chrom_map(chrom_map_str):
    """Parses a chromosome mapping string into a dictionary.

    Args:
        chrom_map_str (str): Mapping string in the format "1=chr1,2=chr2".

    Returns:
        dict: Dictionary mapping VCF chromosome names to FASTA chromosome names.
    """
    if not chrom_map_str:
        return None

    chrom_map = {}
    try:
        for mapping in chrom_map_str.split(","):
            vcf_chrom, fasta_chrom = mapping.split("=")
            chrom_map[vcf_chrom] = fasta_chrom
    except ValueError:
        logger.error("Invalid format for --chrom-map. Use '1=chr1,2=chr2'")
        exit(1)

    return chrom_map

def main():
    """CLI entry point for the vcf_consensus package."""
    parser = argparse.ArgumentParser(description="Generate consensus sequences from VCF and FASTA.")
    parser.add_argument("--vcf", required=True, help="Path to the VCF file.")
    parser.add_argument("--fasta", required=True, help="Path to the FASTA reference genome.")
    parser.add_argument("--length", type=int, required=True, help="Length of consensus sequences.")
    parser.add_argument("--count", type=int, required=True, help="Number of consensus sequences to generate.")
    parser.add_argument("--threshold", type=float, default=0.0, help="Allele frequency threshold for variant inclusion.")
    parser.add_argument("--output", required=True, help="Output FASTA file for consensus sequences.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    parser.add_argument("--chrom-map", type=str, help="Manual chromosome name mapping (e.g., '1=chr1,2=chr2').")
    parser.add_argument("--mode", type=str, choices=["random", "sequential"], default="random",
                        help="Mode of consensus sequence generation: 'random' (default) or 'sequential'.")

    args = parser.parse_args()

    logger.info("Running with parameters:")
    for k, v in vars(args).items():
        logger.info(f"    {k}: {v}")

    chrom_map = parse_chrom_map(args.chrom_map)

    start_time = time.perf_counter()

    try:
        generate_consensus_sequences(
            vcf_path=args.vcf,
            fasta_path=args.fasta,
            length=args.length,
            count=args.count,
            threshold=args.threshold,
            output_path=args.output,
            seed=args.seed,
            chrom_map=chrom_map,
            mode=args.mode
        )
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        exit(1)

    elapsed_time = time.perf_counter() - start_time
    logger.info(f"Finished! Consensus sequences saved in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()