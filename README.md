Проблема связана с тем, что GitHub не распознал Markdown-разметку из-за неправильных отступов или пробелов. Нужно немного поправить README.md, чтобы таблицы и кодовые блоки отображались правильно.

⸻

Исправленный README.md

# vcf_consensus
### Generate Consensus DNA Sequences from VCF and FASTA

`vcf_consensus` is a CLI application written in Python that generates consensus DNA sequences for individual samples based on variations in a VCF file and a reference genome in FASTA format.

Supports both **uncompressed** and **compressed** files (`.vcf.gz`, `.fasta.gz`).

---

## **Installation**
### 1. Clone the repository
```bash
git clone https://github.com/yourusername/vcf_consensus.git
cd vcf_consensus

2. Install in a virtual environment

python -m venv .venv
source .venv/bin/activate  # Linux/Mac
# .venv\Scripts\activate   # Windows

uv pip install -e .



⸻

Usage

Basic execution

python -m vcf_consensus.cli \
    --vcf input.vcf \
    --fasta reference.fasta \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output output_consensus.fasta



⸻

Command-line arguments

Argument	Description
--vcf	Path to the VCF file (supports .vcf.gz).
--fasta	Path to the FASTA reference genome (supports .fasta.gz).
--length	Length of each consensus sequence.
--count	Number of consensus sequences to generate.
--threshold	Allele frequency threshold for variant inclusion (default: 0.0).
--output	Output FASTA file for consensus sequences.
--seed	Random seed for reproducibility (default: None).
--chrom-map	Manually specify chromosome name mapping if VCF and FASTA names do not match (e.g., "1=chr1,2=chr2").



⸻

Examples

1. Basic execution

python -m vcf_consensus.cli \
    --vcf example.vcf \
    --fasta example.fasta \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output consensus.fasta

2. Using .vcf.gz and .fasta.gz

python -m vcf_consensus.cli \
    --vcf example.vcf.gz \
    --fasta example.fasta.gz \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output consensus.fasta

3. Specifying chromosome name mapping

python -m vcf_consensus.cli \
    --vcf example.vcf \
    --fasta example.fasta \
    --length 500 \
    --count 10000 \
    --threshold 0.5 \
    --output consensus.fasta \
    --chrom-map "1=chr1,2=chr2"



⸻

Project Structure

vcf_consensus/
│── vcf_consensus/
│   ├── __init__.py
│   ├── cli.py             # CLI interface
│   ├── consensus.py       # Consensus sequence generation
│   ├── fasta_parser.py    # FASTA parser
│   ├── vcf_parser.py      # VCF parser
│   ├── logger.py          # Logging configuration
│── tests/                 # Unit tests (coming soon)
│── README.md              # Project documentation
│── pyproject.toml         # Package configuration
│── setup.cfg              # Installation configuration (if needed)
│── .gitignore             # Excluded files (virtual environment, logs)



⸻

Development

To modify the project, set up the virtual environment and install the package in editable mode:

uv pip install -e .

Run the test suite before committing changes:

pytest tests/



⸻

License

This project is distributed under the MIT License.

---


If there are any other issues, let me know! 🚀
