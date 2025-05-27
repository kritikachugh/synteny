#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

echo "===== Complete Analysis Pipeline Started ====="
echo "Timestamp: $(date)"
echo "Current working directory: $(pwd)"
echo ""

# --- Configuration ---
SCRIPT_DIR="./code" # Directory where the scripts are located
PYTHON_EXECUTABLE="python3" # Use python3 explicitly, or python if it points to python3

# --- Primary Input Files (Must be present in SCRIPT_DIR) ---
INITIAL_GENE_PRESENCE_ABSENCE_CSV="./data/input/gene_presence_absence.refmt.csv"

# --- Data Directories (Assumed to be subdirectories of SCRIPT_DIR) ---
GFF_DIR="./data/input/gff_files/"
CDS_DIR="./data/input/cds_files/"
GENOME_FASTA_DIR="./data/input/genome_fasta_files/"

# --- Output Directories (Will be created if they don't exist) ---
BLAST_OUTPUT_SUBDIR="./output/blast_final_results/"

# --- Optional: Path to BLAST+ binaries ---
BLAST_PLUS_BIN_PATH="./output/blast_final_results/blast_bin/"  # ADJUST THIS or set to ""

# --- Intermediate and Final Output Filenames ---
# Script 1 outputs
MISSING_GENE_DATA_CSV="./output/missing_gene_data.csv"
UPSTREAM_DOWNSTREAM_CSV="./output/data_upstream_downstream.csv"
# Script 2 outputs
COMPARE_GENES_CSV="./output/compare_genes.csv"
COMPARE_GENES_MATCHES_CSV="./output/compare_genes_matches.csv"
# Script 3 output
COMPARE_GENES_CORE_STATUS_CSV="./output/compare_genes_core_status.csv"
# Script 4 output
ALL_STRAINS_SORTED_WIDE_CSV="./output/all_strains_sorted_with_id.csv"
# Script 5 output
REVERSED_DF_FOR_BLAST_CSV="./output/reversed_df_completeid.csv"
# Script 6 output
SEQUENCE_DATA_CSV="./output/sequence_data.csv"
# Script 7 output
BLAST_NO_HITS_SUMMARY_CSV="./output/blast_no_hits_summary.csv"

# --- Helper function to check for file existence ---
check_file() {
  if [ ! -f "$1" ]; then
    echo "ERROR: Required input file '$1' not found in $(pwd)."
    exit 1
  fi
}
check_dir() {
  if [ ! -d "$1" ]; then
    echo "ERROR: Required input directory '$1' not found relative to $(pwd)."
    echo "Please create it (e.g., '$1') and populate it with necessary files."
    exit 1
  fi
}

# --- Pre-run Checks ---
echo "Performing pre-run checks for initial inputs and data directories..."
check_file "$INITIAL_GENE_PRESENCE_ABSENCE_CSV"
check_dir "$GFF_DIR"
check_dir "$CDS_DIR"
check_dir "$GENOME_FASTA_DIR"
mkdir -p "$BLAST_OUTPUT_SUBDIR"
echo "Pre-run checks completed."
echo ""

# --- Step 1: Extract Upstream/Downstream Gene Context ---
# echo "--- [Step 1/7] Running upstream_downstream_gene_id.py ---"
# $PYTHON_EXECUTABLE "$SCRIPT_DIR/upstream_downstream_gene_id.py" \
#  --presence_absence_csv "$INITIAL_GENE_PRESENCE_ABSENCE_CSV" \
#  --gff_dir "$GFF_DIR" \
#  --missing_gene_out "$MISSING_GENE_DATA_CSV" \
#  --upstream_downstream_out "$UPSTREAM_DOWNSTREAM_CSV"
# echo "Step 1 completed. Outputs: $MISSING_GENE_DATA_CSV, $UPSTREAM_DOWNSTREAM_CSV"
# echo ""

# --- Step 2: Perform Synteny Matching ---
# echo "--- [Step 2/7] Running synteny_match.py ---"
# check_file "$UPSTREAM_DOWNSTREAM_CSV"
# $PYTHON_EXECUTABLE "$SCRIPT_DIR/synteny_match_plus_minus.py" \
#  --input_context_csv "$UPSTREAM_DOWNSTREAM_CSV" \
#  --gff_dir "$GFF_DIR" \
#  --output_compare_csv "$COMPARE_GENES_CSV" \
#  --match_threshold 2
# echo "Step 2 completed. Output: $COMPARE_GENES_CSV"
# echo ""

# # --- Step 3: Add Core Status to Synteny Match Results ---
# echo "--- [Step 3/7] Running add_core_gene_status.py ---"
# check_file "$COMPARE_GENES_CSV"
# check_file "$INITIAL_GENE_PRESENCE_ABSENCE_CSV"
# $PYTHON_EXECUTABLE "$SCRIPT_DIR/add_core_gene_status.py" \
#   --compare_genes_input_csv "$COMPARE_GENES_CSV" \
#   --presence_absence_csv "$INITIAL_GENE_PRESENCE_ABSENCE_CSV" \
#   --output_enriched_csv "$COMPARE_GENES_CORE_STATUS_CSV"
# echo "Step 3 completed. Output: $COMPARE_GENES_CORE_STATUS_CSV"
# echo ""

# --- Step 4: Create Wide Table (all_strains_sorted_with_id.csv) ---
# echo "--- [Step 4/7] Running create_sorted_strains_table.py ---"
# check_file "$UPSTREAM_DOWNSTREAM_CSV"
# check_file "$COMPARE_GENES_CORE_STATUS_CSV"
# $PYTHON_EXECUTABLE "$SCRIPT_DIR/create_sorted_strains_table.py" \
#   --input_upstream_downstream_csv "$UPSTREAM_DOWNSTREAM_CSV" \
#   --input_compare_genes_core_status_csv "$COMPARE_GENES_CORE_STATUS_CSV" \
#   --output_sorted_wide_csv "$ALL_STRAINS_SORTED_WIDE_CSV" \
#   --min_position -5 \
#   --max_position 5
# echo "Step 4 completed. Output: $ALL_STRAINS_SORTED_WIDE_CSV"
# echo ""

# # --- Step 5: Generate/Refine DataFrame for BLAST Input (Conditional Reversal) ---
# echo "--- [Step 5/7] Running reverse_df_logic.py ---"
# check_file "$ALL_STRAINS_SORTED_WIDE_CSV" 
# $PYTHON_EXECUTABLE "$SCRIPT_DIR/reverse_df_logic.py" \
#   --input_sorted_wide_csv "$ALL_STRAINS_SORTED_WIDE_CSV" \
#   --output_reversed_csv "$REVERSED_DF_FOR_BLAST_CSV" \
#   --min_position -5 \
#   --max_position 5
# echo "Step 5 completed. Output: $REVERSED_DF_FOR_BLAST_CSV"
# echo ""

# # --- Step 6: Extract CDS Sequences ---
# echo "--- [Step 6/7] Running cds_sequence.py ---"
# check_file "$INITIAL_GENE_PRESENCE_ABSENCE_CSV"
# $PYTHON_EXECUTABLE "$SCRIPT_DIR/cds_sequence.py" \
#   --presence_absence_csv "$INITIAL_GENE_PRESENCE_ABSENCE_CSV" \
#   --cds_dir "$CDS_DIR" \
#   --output_csv "$SEQUENCE_DATA_CSV"
# echo "Step 6 completed. Output: $SEQUENCE_DATA_CSV"
# echo ""

# # --- Step 7: Perform BLAST Searches ---
# echo "--- [Step 7/7] Running blast_central_genes.py ---"
# check_file "$REVERSED_DF_FOR_BLAST_CSV" 

# BLAST_BIN_ARG_FOR_SCRIPT=""
# if [ -n "$BLAST_PLUS_BIN_PATH" ]; then
#   BLAST_BIN_ARG_FOR_SCRIPT="--blast_bin_path \"$BLAST_PLUS_BIN_PATH\"" 
# fi

# $PYTHON_EXECUTABLE "$SCRIPT_DIR/blast_central_genes.py" \
#   --input_reversed_csv "$REVERSED_DF_FOR_BLAST_CSV" \
#   --cds_files_dir "$CDS_DIR" \
#   --genome_fasta_dir "$GENOME_FASTA_DIR" \
#   --blast_output_dir "$BLAST_OUTPUT_SUBDIR" \
#   --no_hits_output_file "$BLAST_NO_HITS_SUMMARY_CSV" \
#   $BLAST_BIN_ARG_FOR_SCRIPT
# echo "Step 7 completed. BLAST results are in $BLAST_OUTPUT_SUBDIR. No-hits summary: $BLAST_NO_HITS_SUMMARY_CSV"
# echo ""

echo "===== Complete Analysis Pipeline Finished Successfully ====="
echo "Timestamp: $(date)"
