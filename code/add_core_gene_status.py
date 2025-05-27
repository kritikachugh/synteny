# add_core_status_modified.py
import pandas as pd
import argparse

def main():
    # Set up command line arguments
    parser = argparse.ArgumentParser(description="Adds core_status to compare_genes.csv based on gene_presence_absence.refmt.csv")
    parser.add_argument('--compare_genes_input_csv', required=True, 
                      help="Path to the compare_genes.csv file")
    parser.add_argument('--presence_absence_csv', required=True,
                      help="Path to the gene_presence_absence.refmt.csv file")
    parser.add_argument('--output_enriched_csv', required=True,
                      help="Path for output CSV file with added core_status")
    args = parser.parse_args()

    # Read input files
    try:
        compare_genes_df = pd.read_csv(args.compare_genes_input_csv)
        gene_presence_absence_df = pd.read_csv(args.presence_absence_csv)
    except FileNotFoundError as e:
        print(f"Error: Input file not found. {e}")
        return
    except Exception as e:
        print(f"Error reading input CSV files: {e}")
        return

    # Define strain columns
    strain_columns = ["CFT073","ECOR11_PROKKA","ECOR14_PROKKA","ECOR16_PROKKA",
                     "ECOR23_PROKKA","ECOR44_PROKKA","ECOR47_PROKKA","ECOR48_PROKKA",
                     "ECOR49_PROKKA","ECOR4_PROKKA","ECOR50_PROKKA","ECOR51_PROKKA",
                     "ECOR52_PROKKA","ECOR53_PROKKA","ECOR54_PROKKA","ECOR55_PROKKA",
                     "ECOR56_PROKKA","ECOR60_PROKKA","ECOR61_PROKKA","ECOR62_PROKKA",
                     "ECOR64_PROKKA","ECOR66_PROKKA","ECOR6_PROKKA","ECOR8_PROKKA",
                     "KE21","KE24","KE25","KE26","KE4","KE40","KE41","KE46","KE47",
                     "KE48_PROKKA","KE50","KE51","KE54","KE55","KE58","LRPF007",
                     "LRPF62","MG1655","UTI89","W3110"]

    # Create strain dictionary
    strain_dict = {}
    for _, row in gene_presence_absence_df.iterrows():
        for strain in strain_columns:
            value = row[strain]
            if pd.notna(value):
                isolate_count = row['No. isolates']
                if value not in strain_dict:
                    if isolate_count == 44:
                        core_status = "core"
                    elif isolate_count in [42, 43]:
                        core_status = "soft_core"
                    elif isolate_count < 42:
                        core_status = "non_core"
                    else:
                        core_status = "unknown"
                    strain_dict[value] = (row['Gene'], isolate_count, core_status)

    # Add core status to compare_genes DataFrame
    compare_genes_df['core_status'] = 'unknown'
    compare_genes_df['isolate_count'] = None

    for index, row in compare_genes_df.iterrows():
        gene_id = row['ID']
        if gene_id in strain_dict:
            _, isolate_count, core_status = strain_dict[gene_id]
            compare_genes_df.loc[index, 'core_status'] = core_status
            compare_genes_df.loc[index, 'isolate_count'] = isolate_count
        else:
            print(f"Warning: gene_id '{gene_id}' not found in strain_dict")

    # Write output
    try:
        compare_genes_df.to_csv(args.output_enriched_csv, index=False)
        print(f"Successfully wrote {len(compare_genes_df)} rows to {args.output_enriched_csv}")
    except Exception as e:
        print(f"Error writing output file: {e}")

if __name__ == '__main__':
    main()