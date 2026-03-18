# Snippy SNP Pipeline

Automated Snippy variant-calling pipeline: reads alignment → SNP calling → sample renaming → CSV collection.

![License](https://img.shields.io/badge/license-MIT-blue)
![Shell](https://img.shields.io/badge/shell-bash-green)
![Tool](https://img.shields.io/badge/tool-snippy-orange)

## About

This pipeline automates three steps:
1. **Snippy** variant calling on paired-end FASTQ files
2. **Renaming** output folders/files using a CSV mapping table
3. **Collecting** all `.csv` result files into one directory

## Prerequisites

- Bash ≥ 4.0
- [Snippy](https://github.com/tseemann/snippy) installed and in `$PATH`
- Paired-end FASTQ files (`.fastq.gz`)
- Reference genome in `.gbk` format

## Installation
```bash
git clone https://github.com/yourname/snippy-pipeline.git
cd snippy-pipeline
chmod +x run_pipeline.sh
```

## Input Files
```
project/
├── split_output/          # paired FASTQ files
│   ├── sample01_1.fastq.gz
│   ├── sample01_2.fastq.gz
│   └── ...
├── ref.gbk                # reference genome
├── name_mapping.csv       # old_id,new_name
└── run_pipeline.sh
```

`name_mapping.csv` format:
```
old_id,new_name
sample01,StrainA
sample02,StrainB
```

## Usage
```bash
bash run_pipeline.sh
```

Output will be written to:
- `snippy_results/` — per-sample Snippy output (renamed)
- `all_csv_files/`  — all collected `.csv` files

### Key parameters (edit in script)

| Variable | Default | Description |
|---|---|---|
| `IN_DIR` | `split_output` | Input FASTQ directory |
| `REF` | `ref.gbk` | Reference genome |
| `THREADS` | `8` | CPU threads for Snippy |
| `MINCOV` | `10` | Minimum coverage depth |
| `MINFRAC` | `0.9` | Minimum variant allele fraction |

## Output

Each sample in `snippy_results/` contains standard Snippy output files renamed to `<sample_name>.*`. All `.csv` files are copied to `all_csv_files/`.

## License

MIT

# Snippy Variant Filter （All-variants_filtered_without_missense_variants.py）

Filters Snippy CSV output files by removing synonymous and missense variants, retaining only high-impact mutations for downstream analysis.

## What It Does

Reads all `.csv` files from a Snippy output directory, removes rows where the `EFFECT` column contains `synonymous_variant` or `missense_variant`, and saves the filtered results to a new folder.

## Prerequisites

- Python ≥ 3.7
- pandas
```bash
pip install pandas
```

## Input

A folder of Snippy-generated `.csv` files. Each file should contain an `EFFECT` (or `effect`) column.
```
all_csv_files/
├── StrainA.csv
├── StrainB.csv
└── ...
```

## Configuration

Edit the two path variables at the top of the script before running:
```python
INPUT_FOLDER  = "/path/to/all_csv_files/*.csv"
OUTPUT_FOLDER = "/path/to/output_directory"
```

## Usage
```bash
python variant_filter.py
```

Example output:
```
StrainA.csv: saved 12 mutations
StrainB.csv: No residual mutation after filtration.
done! outputs saved in /path/to/output_directory
```

## Filtering Logic

| Excluded effect | Reason |
|---|---|
| `synonymous_variant` | Silent mutation — no amino acid change |
| `missense_variant` | Single amino acid substitution |

All other effects (e.g. `stop_gained`, `frameshift_variant`, `start_lost`) are retained.

## Output

Filtered `.csv` files are saved to `OUTPUT_FOLDER` with the same filenames as the input. Row order and all columns are preserved.

## License

MIT
