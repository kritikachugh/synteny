# add_core_status_modified.py
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Adds core_status to compare_genes.csv based on gene_presence_absence.refmt.csv.")
    parser.add_argument('--compare_genes_input_csv', required=True, help="Path to the compare_genes.csv file (output of synteny matching).")
    parser.add_argument('--presence_absence_csv', required=True, help="Path to the initial gene_presence_absence.refmt.csv file.")
    parser.add_argument('--output_enriched_csv', required=True, help="Path for the output CSV file with added core_status (e.g., compare_genes_core_status.csv).")
    args = parser.parse_args()

    try:
        compare_genes_df = pd.read_csv(args.compare_genes_input_csv)
        gene_presence_absence_df = pd.read_csv(args.presence_absence_csv)
    except FileNotFoundError as e:
        print(f"Error: Input file not found. {e}")
        return
    except Exception as e:
        print(f"Error reading input CSV files: {e}")
        return
    
    if compare_genes_df.empty:
        print(f"Warning: Input file {args.compare_genes_input_csv} is empty. Output will also be empty, but with core_status column.")
        compare_genes_df['core_status'] = 'unknown' 
        compare_genes_df.to_csv(args.output_enriched_csv, index=False)
        print(f"Empty enriched file created at {args.output_enriched_csv}")
        return

    strain_columns_in_pa_df = ["CFT073","ECOR11_PROKKA","ECOR14_PROKKA","ECOR16_PROKKA","ECOR23_PROKKA","ECOR44_PROKKA","ECOR47_PROKKA","ECOR48_PROKKA","ECOR49_PROKKA","ECOR4_PROKKA","ECOR50_PROKKA","ECOR51_PROKKA","ECOR52_PROKKA","ECOR53_PROKKA","ECOR54_PROKKA","ECOR55_PROKKA","ECOR56_PROKKA","ECOR60_PROKKA","ECOR61_PROKKA","ECOR62_PROKKA","ECOR64_PROKKA","ECOR66_PROKKA","ECOR6_PROKKA","ECOR8_PROKKA","KE21","KE24","KE25","KE26","KE4","KE40","KE41","KE46","KE47","KE48_PROKKA","KE50","KE51","KE54","KE55","KE58","LRPF007","LRPF62","MG1655","UTI89","W3110"]

    roary_id_to_core_status_map = {}
    for _, row in gene_presence_absence_df.iterrows():
        isolate_count = row.get('No. isolates', 0)
        core_status = "unknown"
        if isolate_count == 44: core_status = "core"
        elif isolate_count == 42 or isolate_count == 43: core_status = "soft_core"
        elif isolate_count < 42: core_status = "non_core"
        
        for strain_col_name in strain_columns_in_pa_df:
            if strain_col_name in row and pd.notna(row[strain_col_name]):
                gene_ids_in_cell = str(row[strain_col_name]).split(';')
                for roary_gene_id in gene_ids_in_cell:
                    roary_gene_id = roary_gene_id.strip()
                    if roary_gene_id and roary_gene_id not in roary_id_to_core_status_map:
                        roary_id_to_core_status_map[roary_gene_id] = core_status
                        
    if 'ID' not in compare_genes_df.columns:
        print(f"Error: 'ID' column not found in {args.compare_genes_input_csv}. Cannot map core status.")
        compare_genes_df.to_csv(args.output_enriched_csv, index=False) 
        return

    compare_genes_df['core_status'] = compare_genes_df['ID'].map(roary_id_to_core_status_map).fillna('unknown')
    
    try:
        compare_genes_df.to_csv(args.output_enriched_csv, index=False)
        print(f"Enriched compare_genes data with core status saved to {args.output_enriched_csv}")
    except Exception as e:
        print(f"Error writing output CSV to {args.output_enriched_csv}: {e}")

if __name__ == '__main__':
    main()