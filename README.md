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
bash run_snippy_pipeline.sh
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
| `BASEQUAL` | `20` | Minimum base quality |
| `MAPQUAL` | `30` | Minimum mapping qualit |

## Output

Each sample in `snippy_results/` contains standard Snippy output files renamed to `<sample_name>.*`. All `.csv` files are copied to `all_csv_files/`.


# GenBank CDS Extractor

Parses a GenBank (`.gbk`) reference file and extracts all CDS features into a structured CSV table.

## What It Does

Iterates over every `CDS` feature in a GenBank file and writes the following fields to a CSV:

| Column | Description |
|---|---|
| `Start` | 1-based genomic start position |
| `End` | Genomic end position |
| `Strand` | `+` or `-` |
| `Gene_Symbol` | Gene name (e.g. `mfsA`), `-` if absent |
| `Locus_Tag` | Locus tag (e.g. `BPSL0001`), `-` if absent |
| `Protein_ID` | NCBI protein accession (e.g. `WP_012345678.1`), `-` if absent |
| `Product_Function` | Annotated product description |

## Configuration

Edit the two path variables at the top of the script:
```python
GBK_FILE   = "/path/to/ref.gbk"
OUTPUT_CSV = "/path/to/ref_gene_functions_with_protein_id.csv"
```

## Usage
```bash
python extract_info.py
```

Example terminal output:
```
Done! Extracted 5,832 genes. Output saved to ref_gene_functions_with_protein_id.csv
```

## Output

A CSV file with one row per CDS feature:
```
Start,End,Strand,Gene_Symbol,Locus_Tag,Protein_ID,Product_Function
1423,2301,+,mfsA,BPSL0001,WP_012345678.1,MFS transporter
2450,3102,-,-,BPSL0002,WP_012345679.1,hypothetical protein
```

# Snippy Variant Filter 

Filters Snippy CSV output files by removing synonymous and missense variants, retaining only high-impact mutations for downstream analysis.

## What It Does

Reads all `.csv` files from a Snippy output directory, removes rows where the `EFFECT` column contains `synonymous_variant` or `missense_variant`, and saves the filtered results to a new folder.

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
python All-variants_filtered_without_missense_variants.py
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


# MFS Superfamily Mutation Matrix

Builds a gene × sample mutation count matrix for a target gene list (e.g. MFS superfamily locus tags) from filtered Snippy variant CSV files.

## What It Does

1. Reads a list of target locus tags from `MFS_mutation.csv`
2. Scans all filtered variant `.csv` files for hits matching those tags
3. Outputs a pivot matrix — rows = locus tags, columns = sample names, values = mutation counts


## Project Structure
```
project_1b/
├── All_variants_filtered_without_missense_variants/
│   ├── StrainA.csv
│   ├── StrainB.csv
│   └── ...
└── MFS_superfamily/
    ├── MFS_mutation.csv        ← target locus tag list (input)
    └── MFS_Mutation_Matrix.csv ← output matrix
```

## Input Files

**`MFS_mutation.csv`** — one locus tag per row, no header required (auto-detected)

**Variant CSV files** — output from the Snippy Variant Filter step, must contain a `LOCUS_TAG` (or `locus_tag`) column.

## Configuration

Edit `BASE_DIR` at the top of the script to match your local path:
```python
BASE_DIR = "./project_1b"
```

All other paths are derived automatically.

## Usage
```bash
python Target_Gene_mutation_Matrix.py
```

Example output:
```
Loaded 38 locus tags from first column.
Found 45 sample files, start screening and counting...
All done! The mutation matrix for target genes has been generated.
Target gene (row):    38
Sample size (column): 45
```

## Output

`MFS_Mutation_Matrix.csv` — a gene × sample count matrix:

|  | StrainA | StrainB | StrainC |
|---|---|---|---|
| PMI_RS0001 | 2 | 0 | 1 |
| PMI_RS0042 | 0 | 3 | 0 |
| PMI_RS1234 | 1 | 1 | 2 |

- Rows are sorted alphabetically by locus tag
- Columns are sorted alphabetically by sample name
- Locus tags with zero mutations across all samples are included (value = `0`)
- `.csv` extension is stripped from column headers

## Pipeline Position
```
Snippy CSV output
      ↓
variant_filter.py       ← removes synonymous/missense variants
      ↓
mutation_matrix.py      ← this script
      ↓
MFS_Mutation_Matrix.csv
```


# MFS Mutation Heatmap

Visualises a gene × sample mutation count matrix as a paginated, group-annotated heatmap, with samples ordered by CHD MIC level and colour-coded by MIC group.

## What It Does

1. Loads a mutation count matrix (`MFS_Mutation_Matrix.csv`)
2. Orders samples according to a strain order file (high → low CHD MIC)
3. Annotates row labels with gene descriptions from `MFS_mutation.csv`
4. Renders a paginated PDF heatmap with MIC group colour bars above the columns

## Example Output

| Feature | Detail |
|---|---|
| Rows | Locus tags + gene descriptions |
| Columns | Sample names, ordered by MIC |
| Cell values | Mutation counts (capped at `CUSTOM_VMAX`) |
| Group bar | Colour-coded MIC range above columns |
| Format | PDF, one file per page |

##  Structure
```
_1b/
├── strain_order_CHD_MICs_hightolow.csv   ← sample order + MIC group labels
├── MFS_superfamily/
│   ├── MFS_mutation.csv                  ← locus tags + gene descriptions
│   ├── MFS_Mutation_Matrix.csv           ← input matrix (from mutation_matrix.py)
│   └── Heatmap/
│       └── Heatmap_Grouped_by_CHD_MIC_grouped.pdf  ← output
```

## Input Files

**`strain_order_CHD_MICs_hightolow.csv`** — defines column order and MIC grouping:
```
MIC range, Strains
High MIC,  StrainA
High MIC,  StrainB
Low MIC,   StrainC
```

- `MIC range` column: group label (colour-coded; any name containing `high` → red, `low` → green, other → yellow)
- `Strains` column: sample names matching matrix column headers
- Empty cells in `MIC range` are forward-filled automatically

**`MFS_mutation.csv`** — locus tags with optional gene descriptions (second column):
```
PMI_RS0001, MFS transporter
PMI_RS0042, TetR family regulator
```

## Configuration

Edit the constants at the top of the script:
```python
BASE_DIR       = "./project_1b"
GENES_PER_PAGE = 100     # rows per PDF page
CUSTOM_VMAX    = 3       # colour scale ceiling
CELL_HEIGHT    = 0.22    # inches per row
CELL_WIDTH     = 0.22    # inches per column
```

Colours use a custom 4-stop blue gradient (`#FAFBFD` → `#2B3F5E`).

## Usage

**Basic:**
```bash
python heatmap_with_annotation.py
```

**With CLI overrides:**
```bash
python heatmap_with_annotation.py \
  --cell-height 0.3 \
  --cell-width 0.25 \
  --genes-per-page 50
```

### CLI Arguments

| Argument | Default | Description |
|---|---|---|
| `--cell-height` | `0.22` | Row height in inches |
| `--cell-width` | `0.22` | Column width in inches |
| `--genes-per-page` | `100` | Genes (rows) per PDF page |

## Pipeline Position
```
MFS_Mutation_Matrix.csv
         ↓
    heatmap.py          ← this script
         ↓
Heatmap_Grouped_by_CHD_MIC_grouped.pdf
```

# PCoA Analysis of MFS Regulator Mutations

Performs Principal Coordinates Analysis (PCoA) on a binary mutation matrix using Jaccard distance, coloured by CHD MIC group, with PERMANOVA significance testing.

## What It Does

1. Loads a gene × sample mutation matrix and CHD MIC data
2. Computes pairwise **Jaccard distances** between samples
3. Runs **PCoA** (classical MDS) to project samples into 2D space
4. Tests group separation with **PERMANOVA** (999 permutations)
5. Plots scatter with **95% confidence ellipses** per MIC group and strain labels for High-MIC samples
6. Saves output as a PDF


| Feature | Detail |
|---|---|
| Distance metric | Jaccard |
| Axes | PCoA 1 & 2 with % variance explained |
| Colours | Red = High MIC, Yellow = Intermediate, Green = Low |
| Ellipses | 1.8 SD confidence ellipse per group (≥ 3 samples) |
| Labels | Strain names shown for High-MIC group only |
| Statistics | PERMANOVA pseudo-F and p-value in plot corner |


## Project Structure
```
project_1b/
├── MIC_data.csv                              ← CHD MIC values per strain
└── MFS_selected_genes/
    ├── MFS_regulator_mutation_Matrix.csv     ← input mutation matrix
    └── PCoA_analysis/
        └── PCoA_MFS_Regulator.pdf            ← output plot
```

## Input Files

**`MFS_regulator_mutation_Matrix.csv`** — gene × sample or sample × gene mutation count matrix (orientation is detected automatically):
```
         StrainA  StrainB  StrainC
PMI_RS0001       1        0        2
PMI_RS0042       0        1        0
```

**`MIC_data.csv`** — CHD MIC values per strain (column containing `chd` is auto-detected):
```
Strains,  CHD_MIC
StrainA,  8
StrainB,  4
StrainC,  2
```

- Index column is auto-detected from: `strains`, `strain`, or `sample`
- `.csv` suffixes in strain names are stripped automatically

## MIC Group Thresholds

Edit at the top of the script:
```python
THRESHOLD_HIGH = 7.0   # MIC ≥ 7  → High (red)
THRESHOLD_LOW  = 4.0   # MIC ≤ 4  → Low  (green)
                       # between  → Intermediate (yellow)
```

## Configuration
```python
BASE_DIR      = "./project_1b"
OUTPUT_PLOT   = "...PCoA_MFS_Regulator.pdf"
```

## Usage
```bash
python PCoA_analysis.py
```

Example terminal output:
```
All done! Saved in /path/to/PCoA_analysis/PCoA_MFS_Regulator.pdf
```

## Method Notes

| Step | Detail |
|---|---|
| Distance | Jaccard (binary mutation presence/absence) |
| PCoA | Classical MDS via eigen-decomposition of centred distance matrix |
| Variance | Calculated from positive eigenvalues only |
| Jitter | Small Gaussian noise added to avoid overplotting |
| PERMANOVA | 999 permutations, random seed fixed at 42 for reproducibility |
| Ellipses | Pearson correlation ellipse scaled to 1.8 SD; requires ≥ 3 samples per group |

## Pipeline Position
```
MFS_regulator_mutation_Matrix.csv
          +
      MIC_data.csv
          ↓
    pcoa_analysis.py        ← this script
          ↓
PCoA_MFS_Regulator.pdf
```

# Elastic Net Logistic Regression — CHD Tolerance Feature Selection

Trains Elastic Net logistic regression models across a range of L1/L2 penalty ratios to identify MFS gene mutations associated with high CHD (chlorhexidine digluconate) MIC, with cross-validation accuracy and feature sparsity comparison.

## What It Does

1. Loads a gene × sample mutation matrix and CHD MIC data
2. Filters out genes with fewer than `MIN_MUTATION_COUNT` mutations across all samples
3. Binarises samples into **High** vs **Low** MIC groups
4. Fits **Elastic Net logistic regression** at 7 L1 ratio values (pure Ridge → pure LASSO)
5. Evaluates each model with **5×5 RepeatedStratifiedKFold** cross-validation
6. Outputs per-ratio feature importance bar plots and an L1 ratio comparison panel

## Output Files

| File | Description |
|---|---|
| `ElasticNet_Feature_Importance_CHD_l1ratio_*.pdf` | Top 10 features per L1 ratio |
| `ElasticNet_L1_ratio_Comparison.pdf` | CV accuracy + active feature count across all ratios |


## Project Structure
```
project_1b/
├── MIC_data.csv
└── MFS_selected_genes/
    ├── MFS_selected_genes_mutation_Matrix.csv   ← input mutation matrix
    ├── MFS_selected_genes.csv                   ← locus tag annotations
    └── Elastic_Net_Outcomes_final/              ← output plots
```

## Input Files

**`MFS_selected_genes_mutation_Matrix.csv`** — gene × sample or sample × gene mutation count matrix (orientation auto-detected):
```
         StrainA  StrainB  StrainC
PMI_RS0001       1        0        2
PMI_RS0042       0        1        0
```

**`MIC_data.csv`** — CHD MIC values per strain (column containing `chd` is auto-detected):
```
Strains,  CHD_MIC
StrainA,  8
StrainB,  4
StrainC,  2
```

**`MFS_selected_genes.csv`** — locus tag to gene description mapping (no header):
```
PMI_RS0001, MFS transporter
PMI_RS0042, TetR family regulator
```

## Configuration
```python
THRESHOLD_HIGH     = 7.0   # MIC ≥ 7  → High (positive class)
THRESHOLD_LOW      = 4.0   # MIC ≤ 4  → Low  (negative class)
                            # between  → excluded from analysis
MIN_MUTATION_COUNT = 3      # minimum mutations across all samples to retain a gene
L1_ratioS = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
```

## Usage
```bash
python elasticNet_L1_Ratio_Comparison.py
```

Example terminal output:
```
Loaded 38 annotations.
Best L1_ratio = 0.5 (ElasticNet) with CV Accuracy = 0.81
```

## Method Notes

| Step | Detail |
|---|---|
| Scaling | `StandardScaler` (zero mean, unit variance) |
| Solver | `saga` (supports all penalty combinations) |
| Class weight | `balanced` (corrects for High/Low imbalance) |
| Cross-validation | `RepeatedStratifiedKFold` — 5 splits × 5 repeats |
| Error bars | SEM = CV std / √n |
| Regularisation strength | Fixed at `C = 0.5` across all ratios |
| Feature threshold | Coefficients < 1×10⁻⁵ treated as zero (inactive) |
| Reproducibility | `random_state = 42` throughout |

### L1 Ratio Guide

| L1 ratio | Penalty | Behaviour |
|---|---|---|
| `0.0` | Pure Ridge (L2) | Shrinks all coefficients, retains all features |
| `0.1–0.9` | Elastic Net | Balances sparsity and grouping of correlated features |
| `1.0` | Pure LASSO (L1) | Maximum sparsity — drives many coefficients to exactly zero |

### Interpreting Coefficients

- **Positive coefficient** → mutation associated with **resistance** (High MIC)
- **Negative coefficient** → mutation associated with **susceptibility** (Low MIC)

## Pipeline Position
```
MFS_selected_genes_mutation_Matrix.csv
              +
          MIC_data.csv
              ↓
    elastic_net_analysis.py       ← this script
              ↓
  ElasticNet_Feature_Importance_CHD_l1ratio_*.pdf
  ElasticNet_L1_ratio_Comparison.pdf
```

## License

MIT
