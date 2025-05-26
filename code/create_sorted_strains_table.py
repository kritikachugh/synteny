# create_sorted_strains_table.py
import pandas as pd
import argparse

def process_group_to_wide_format(df_group, group_central_gene_name, min_pos=-5, max_pos=5):
    reshaped_rows_list = []
    for strain_val, strain_specific_data in df_group.groupby('strain'):
        positions_data_map = {}
        for _, row_item in strain_specific_data.iterrows():
            position_val = row_item.get('position')
            gene_tuple_info = (
                row_item.get('gene'),
                row_item.get('ID'),
                row_item.get('strand'),
                row_item.get('core_status'),
                row_item.get('start'),
                row_item.get('end')
            )
            if pd.notna(position_val):
                 positions_data_map[int(position_val)] = gene_tuple_info
        
        output_row_dict = {'central_gene': group_central_gene_name, 'strain': strain_val}
        for i in range(min_pos, max_pos + 1):
            col_name_pos = f"position_{i}"
            output_row_dict[col_name_pos] = positions_data_map.get(i) 
        reshaped_rows_list.append(output_row_dict)
    return pd.DataFrame(reshaped_rows_list)

def main():
    parser = argparse.ArgumentParser(description="Reshapes and combines gene context data into a single sorted wide-format CSV.")
    parser.add_argument('--input_upstream_downstream_csv', required=True, help="Path to data_upstream_downstream.csv.")
    parser.add_argument('--input_compare_genes_core_status_csv', required=True, help="Path to compare_genes_core_status.csv.")
    parser.add_argument('--output_sorted_wide_csv', required=True, help="Output path for the sorted wide format CSV (e.g., all_strains_sorted_with_id.csv).")
    parser.add_argument('--min_position', type=int, default=-5, help="Minimum relative position (default: -5).")
    parser.add_argument('--max_position', type=int, default=5, help="Maximum relative position (default: 5).")
    args = parser.parse_args()

    try:
        df_upstream_downstream = pd.read_csv(args.input_upstream_downstream_csv)
        df_compare_genes = pd.read_csv(args.input_compare_genes_core_status_csv) 
    except FileNotFoundError as e:
        print(f"Error: Input file not found. {e}")
        return
    except Exception as e:
        print(f"Error reading input CSV files: {e}")
        return

    required_cols = ['central_gene', 'strain', 'position', 'gene', 'ID', 'strand', 'core_status', 'start', 'end']
    for df, df_name in [(df_upstream_downstream, args.input_upstream_downstream_csv), 
                        (df_compare_genes, args.input_compare_genes_core_status_csv)]:
        if df.empty: # Allow empty input DFs, will result in empty section in output
            print(f"Warning: Input DataFrame from '{df_name}' is empty.")
            continue
        for col in required_cols:
            if col not in df.columns:
                print(f"Error: Required column '{col}' missing in DataFrame from '{df_name}'. Cannot proceed.")
                return
    
    all_strains_df_wide_list = []
    if not df_upstream_downstream.empty:
        print("Processing 'data_upstream_downstream.csv' into wide format...")
        all_strains_df_wide_list = df_upstream_downstream.groupby('central_gene', group_keys=False).apply(
            lambda x: process_group_to_wide_format(x, x.name, args.min_position, args.max_position)
        ).reset_index(drop=True)
    else:
        all_strains_df_wide_list = pd.DataFrame()


    missing_strains_df_wide_list = []
    if not df_compare_genes.empty:
        print("Processing 'compare_genes_core_status.csv' into wide format...")
        missing_strains_df_wide_list = df_compare_genes.groupby('central_gene', group_keys=False).apply(
            lambda x: process_group_to_wide_format(x, x.name, args.min_position, args.max_position)
        ).reset_index(drop=True)
    else:
        missing_strains_df_wide_list = pd.DataFrame()

    # Add priority for sorting
    if not missing_strains_df_wide_list.empty:
        missing_strains_df_wide_list["priority"] = 0 
    if not all_strains_df_wide_list.empty:
        all_strains_df_wide_list["priority"] = 1  
    
    # Define columns for empty dataframes if one input was empty
    position_cols = [f"position_{i}" for i in range(args.min_position, args.max_position + 1)]
    base_cols = ['central_gene', 'strain'] + position_cols
    
    if missing_strains_df_wide_list.empty and not isinstance(missing_strains_df_wide_list, pd.DataFrame):
        missing_strains_df_wide_list = pd.DataFrame(columns=base_cols + ['priority'])
    if all_strains_df_wide_list.empty and not isinstance(all_strains_df_wide_list, pd.DataFrame):
         all_strains_df_wide_list = pd.DataFrame(columns=base_cols + ['priority'])


    final_sorted_df = pd.DataFrame()
    if not missing_strains_df_wide_list.empty or not all_strains_df_wide_list.empty:
        combined_df_for_sort = pd.concat([missing_strains_df_wide_list, all_strains_df_wide_list], ignore_index=True)
        if not combined_df_for_sort.empty and 'central_gene' in combined_df_for_sort.columns:
            final_sorted_df = combined_df_for_sort.sort_values(by=['central_gene', 'priority'])
            if 'priority' in final_sorted_df.columns:
                 final_sorted_df = final_sorted_df.drop(columns=['priority'])
        else: # Handle if combined is still empty or missing sort columns
            final_sorted_df = combined_df_for_sort.drop(columns=['priority'], errors='ignore')

    try:
        final_sorted_df.to_csv(args.output_sorted_wide_csv, index=False)
        print(f"Combined, sorted wide format data ({len(final_sorted_df)} rows) saved to {args.output_sorted_wide_csv}")
    except Exception as e:
        print(f"Error writing output CSV {args.output_sorted_wide_csv}: {e}")

if __name__ == '__main__':
    main()