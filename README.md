# vcf_consensus
### Generate Consensus DNA Sequences from VCF and FASTA

`vcf_consensus` is a CLI application written in Python that generates consensus DNA sequences for individual samples based on variations in a VCF file and a reference genome in FASTA format.

Supports both **uncompressed** and **compressed** files (`.vcf.gz`, `.fasta.gz`).

---

## **Installation**

### **!Make sure that you have uv installed in your system, current version uses it as a package manager!**

```bash
pip install uv 
```
(in case you don't have one)

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/vcf_consensus.git
cd vcf_consensus
```

### 2. Install in a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate  # Linux/Mac
# .venv\Scripts\activate   # Windows

uv pip install -e .
```

## **Usage**

### Basic execution
```bash
python -m vcf_consensus.cli \
    --vcf input.vcf \
    --fasta reference.fasta \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output output_consensus.fasta
```

### Command-line options
```
--vcf	        Path to the VCF file (supports .vcf.gz).
--fasta	        Path to the FASTA reference genome (supports .fasta.gz).
--length	    Length of each consensus sequence.
--count	        Number of consensus sequences to generate.
--threshold     Probability of using an alternative allele when present in the VCF (default: 0.0). 
                A value of 0.5 means that 50% of heterozygous positions will use an alternative allele.
--output	    Output FASTA file for consensus sequences.
--seed	        Random seed for reproducibility (default: None).
--chrom-map     (Optinal) Manual chromosome name mapping if VCF and FASTA names differ (e.g., "1=chr1,2=chr2").
                Use this only if chromosomes in the VCF do not match those in the FASTA.
--mode          Method for selecting consensus start positions (default: random). 
                Options:
                  - "random": Randomly selects start positions.
                  - "sequential": Selects positions sequentially with 50% overlap.
```


## Examples

### 1. Basic execution

```bash
python -m vcf_consensus.cli \
    --vcf example.vcf \
    --fasta example.fasta \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output consensus.fasta
```

### 2. Using .vcf.gz and .fasta.gz

```bash
python -m vcf_consensus.cli \
    --vcf example.vcf.gz \
    --fasta example.fasta.gz \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output consensus.fasta
```

### 3. Specifying chromosome name mapping

```bash
python -m vcf_consensus.cli \
    --vcf example.vcf \
    --fasta example.fasta \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output consensus.fasta \
    --chrom-map "1=chr1,2=chr2"
```

## Performance Considerations

- Large FASTA files (>3 GB) are fully loaded into memory. Might change later
- Large VCF files (>10 million variants) may cause high memory usage.
- If you experience slow performance, reduce --count or increase --length to reduce the number of consensus sequences.

License

This project is distributed under the MIT License.
