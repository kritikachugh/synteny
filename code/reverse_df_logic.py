# reverse_df_logic.py
import pandas as pd
import argparse
import ast 

def get_gene_id_from_pos_tuple_string(cell_value_str):
    if pd.isna(cell_value_str):
        return None
    try:
        actual_tuple = ast.literal_eval(str(cell_value_str))
        if isinstance(actual_tuple, tuple) and len(actual_tuple) > 0:
            return str(actual_tuple[0]).strip() 
    except (ValueError, SyntaxError):
        return str(cell_value_str).strip()
    return str(cell_value_str).strip()

def reverse_columns_in_row(row_series, min_pos=-5, max_pos=5): # Made range configurable
    cols_to_reverse_numeric_range = range(min_pos, max_pos + 1) 

    cols_to_reverse = []
    for i in cols_to_reverse_numeric_range:
        if i == 0: continue 
        col_name = f"position_{i}"
        if col_name in row_series.index:
            cols_to_reverse.append(col_name)
            
    if not cols_to_reverse:
        return row_series

    cols_to_reverse.sort(key=lambda x: int(x.split('_')[1]))
    
    new_row_data = row_series.copy()
    current_values_in_order = [row_series[col] for col in cols_to_reverse]
    reversed_values_for_assignment = current_values_in_order[::-1]

    for i, col_name_ordered in enumerate(cols_to_reverse):
        new_row_data[col_name_ordered] = reversed_values_for_assignment[i]
        
    return new_row_data

def process_dataframe_for_conditional_reversal(df_input, min_pos=-5, max_pos=5):
    df_processed = df_input.copy()
    
    if 'central_gene' not in df_processed.columns or 'position_0' not in df_processed.columns:
        print("Error: 'central_gene' or 'position_0' column missing. Cannot process for reversal.")
        return df_input

    for central_gene_name, group_df in df_processed.groupby('central_gene'):
        empty_pos0_rows = group_df[group_df['position_0'].isnull()] 
        non_empty_pos0_rows = group_df[group_df['position_0'].notnull()]

        if empty_pos0_rows.empty or non_empty_pos0_rows.empty:
            continue
        
        empty_pos0_template_row = empty_pos0_rows.iloc[0]
        empty_pos0_positive_set = set()
        empty_pos0_negative_set = set()

        for i in range(min_pos, max_pos + 1): 
            if i == 0: continue 
            key_col = f"position_{i}"
            if key_col in empty_pos0_template_row.index:
                cell_val_str = empty_pos0_template_row.get(key_col)
                if pd.notna(cell_val_str):
                    gene_id = get_gene_id_from_pos_tuple_string(cell_val_str)
                    if gene_id:
                        if i > 0: empty_pos0_positive_set.add(gene_id)
                        elif i < 0: empty_pos0_negative_set.add(gene_id)
        
        for row_idx, non_empty_pos0_row_original in non_empty_pos0_rows.iterrows():
            non_empty_pos0_positive_set = set()
            non_empty_pos0_negative_set = set()
            for i in range(min_pos, max_pos + 1):
                if i == 0: continue
                key_col = f"position_{i}"
                if key_col in non_empty_pos0_row_original.index:
                    cell_val_str = non_empty_pos0_row_original.get(key_col)
                    if pd.notna(cell_val_str):
                        gene_id = get_gene_id_from_pos_tuple_string(cell_val_str)
                        if gene_id:
                            if i > 0: non_empty_pos0_positive_set.add(gene_id)
                            elif i < 0: non_empty_pos0_negative_set.add(gene_id)
            
            original_match_count = len(empty_pos0_positive_set.intersection(non_empty_pos0_positive_set)) + \
                                   len(empty_pos0_negative_set.intersection(non_empty_pos0_negative_set))
            reverse_match_count = len(empty_pos0_positive_set.intersection(non_empty_pos0_negative_set)) + \
                                  len(empty_pos0_negative_set.intersection(non_empty_pos0_positive_set))

            if reverse_match_count > original_match_count:
                print(f"Reversing position columns for central_gene '{central_gene_name}', original row index {row_idx} (strain: {non_empty_pos0_row_original.get('strain', 'N/A')}).")
                df_processed.loc[row_idx] = reverse_columns_in_row(non_empty_pos0_row_original.copy(), min_pos, max_pos)
    return df_processed

def main():
    parser = argparse.ArgumentParser(description="Conditionally reverses positional gene columns in a CSV.")
    parser.add_argument('--input_sorted_wide_csv', required=True, help="Input CSV (e.g., all_strains_sorted_with_id.csv).")
    parser.add_argument('--output_reversed_csv', required=True, help="Output CSV (e.g., reversed_df_completeid.csv).")
    parser.add_argument('--min_position', type=int, default=-5, help="Minimum position for reversal range (e.g. -5).")
    parser.add_argument('--max_position', type=int, default=5, help="Maximum position for reversal range (e.g. 5).")
    args = parser.parse_args()

    try:
        df_input_wide = pd.read_csv(args.input_sorted_wide_csv)
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_sorted_wide_csv}' not found.")
        return
    except Exception as e:
        print(f"Error reading input CSV {args.input_sorted_wide_csv}: {e}")
        return

    if df_input_wide.empty:
        print(f"Input file {args.input_sorted_wide_csv} is empty. Writing an empty output file.")
        df_input_wide.to_csv(args.output_reversed_csv, index=False)
        return

    print(f"Processing {args.input_sorted_wide_csv} for conditional reversal...")
    df_conditionally_reversed = process_dataframe_for_conditional_reversal(df_input_wide, args.min_position, args.max_position)
    
    try:
        df_conditionally_reversed.to_csv(args.output_reversed_csv, index=False)
        print(f"Conditionally reversed data saved to {args.output_reversed_csv}")
    except Exception as e:
        print(f"Error writing output CSV {args.output_reversed_csv}: {e}")

if __name__ == '__main__':
    main()