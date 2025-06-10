# **Gene Context Analysis Pipeline**

This repository contains Python scripts for analyzing gene context, synteny, and core gene status in bacterial genomes. The pipeline integrates data from Roary, BLAST, and GFF files to provide insights into gene presence/absence, synteny, and core/non-core gene classification.

---

## **Table of Contents**
- Overview
- Features
- Requirements
- Installation
- Usage
- Scripts
  - 1. `synteny_match_plus_minus.py`
  - 2. `upstream_downstream_gene_id.py`
  - 3. `blast_central_genes.py`
  - 4. `add_core_gene_status.py`
- Input Files
- Output Files
- Troubleshooting
- License

---

## **Overview**
This pipeline is designed to:
1. Extract upstream and downstream gene context for central genes.
2. Perform BLAST analysis to identify missing genes in specific strains.
3. Add core/non-core gene status to the results based on Roary's gene presence/absence data.

---

## **Features**
- Parse GFF files to extract gene features.
- Identify syntenic regions for missing genes.
- Perform BLAST-based searches for missing genes.
- Classify genes as core, soft-core, or non-core based on isolate counts.
- Generate enriched CSV outputs for downstream analysis.

---

## **Requirements**
- Python 3.8 or higher
- Required Python libraries:
  - `pandas`
  - `argparse`
  - `tqdm`
- BLAST+ tools (`makeblastdb`, `blastn`)

---

## **Installation**
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/gene-context-analysis.git
   cd gene-context-analysis
   ```

2. Install required Python libraries:
   ```bash
   pip install -r requirements.txt
   ```

3. Ensure BLAST+ tools are installed and accessible in your system's PATH.

---

## **Usage**
Each script in the repository is designed for a specific task. Below are the usage instructions for each script.

---

## **Scripts**

### **1. `synteny_match_plus_minus.py`**
**Description**: Identifies syntenic regions for missing genes among strains.

**Usage**:
```bash
python synteny_match_plus_minus.py --input_context_csv <input_csv> --gff_dir <gff_directory> --output_compare_csv <output_csv> --output_matches_csv <matches_csv>
```

**Arguments**:
- `--input_context_csv`: Path to the upstream/downstream context CSV.
- `--gff_dir`: Directory containing GFF files.
- `--output_compare_csv`: Output CSV for detailed syntenic regions.
- `--output_matches_csv`: Output CSV for summary of matches.

---

### **2. `upstream_downstream_gene_id.py`**
**Description**: Extracts upstream and downstream gene IDs for central genes.

**Usage**:
```bash
python upstream_downstream_gene_id.py --presence_absence_csv <presence_absence_csv> --gff_dir <gff_directory> --output_csv <output_csv>
```

**Arguments**:
- `--presence_absence_csv`: Path to the gene presence/absence file.
- `--gff_dir`: Directory containing GFF files.
- `--output_csv`: Output CSV file.

---

### **3. `blast_central_genes.py`**
**Description**: Performs BLAST analysis to identify missing genes in specific strains.

**Usage**:
```bash
python blast_central_genes.py --input_reversed_csv <input_csv> --genome_fasta_dir <fasta_directory> --output_dir <output_directory>
```

**Arguments**:
- `--input_reversed_csv`: Path to the reversed input CSV.
- `--genome_fasta_dir`: Directory containing genome FASTA files.
- `--output_dir`: Directory for BLAST results.

---

### **4. add_core_gene_status.py**
**Description**: Adds core/non-core gene status to the `compare_genes.csv` file.

**Usage**:
```bash
python add_core_gene_status.py --compare_genes_input_csv <compare_genes_csv> --presence_absence_csv <presence_absence_csv> --output_enriched_csv <output_csv>
```

**Arguments**:
- `--compare_genes_input_csv`: Path to the `compare_genes.csv` file.
- `--presence_absence_csv`: Path to the `gene_presence_absence.refmt.csv` file.
- `--output_enriched_csv`: Path for the enriched output CSV.

---

## **Input Files**
- **GFF Files**: Genome annotation files for each strain.
- **Gene Presence/Absence File**: Roary's `gene_presence_absence.refmt.csv`.
- **Genome FASTA Files**: FASTA files for each strain.

---

## **Output Files**
- **`compare_genes.csv`**: Detailed syntenic regions for central genes.
- **compare_genes_matches.csv**: Summary of matches for syntenic regions.
- **BLAST Results**: Output files for BLAST analysis.
- **Enriched CSV**: CSV file with added core/non-core gene status.

---

## **Troubleshooting**
1. **UnicodeEncodeError**:
   - Ensure all file operations use `encoding='utf-8'`.

2. **FileNotFoundError**:
   - Verify that input files and directories exist and are correctly specified.

3. **BLAST Errors**:
   - Ensure BLAST+ tools are installed and accessible in your system's PATH.

4. **Missing Genes in Output**:
   - Check if the gene IDs in the input files match those in the GFF files.

---
