# synteny_match_plus_minus.py
import pandas as pd
import os
import csv
import argparse
import re

def parse_gff_to_df(gff_file_path):
    """Parse GFF file into DataFrame using only existing IDs."""
    features = []
    try:
        with open(gff_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    seqid, source, ftype, start, end, score, strand, phase, attributes = fields[:9]
                    attr_dict = {}
                    for attr_pair in attributes.split(';'):
                        if '=' in attr_pair:
                            key, value = attr_pair.split('=', 1)
                            attr_dict[key] = value.strip().strip('"')
                    
                    # Use only existing IDs or locus_tags
                    feature_id = None
                    if 'locus_tag' in attr_dict:
                        feature_id = attr_dict['locus_tag']
                    elif 'ID' in attr_dict:
                        feature_id = attr_dict['ID']
                    
                    feature = {
                        'seqid': seqid,
                        'source': source,
                        'type': ftype,
                        'start': int(start),
                        'end': int(end),
                        'score': score,
                        'strand': strand,
                        'phase': phase,
                        'ID': feature_id,
                        **attr_dict
                    }
                    features.append(feature)
                    
    except FileNotFoundError:
        print(f"Error: GFF file not found: {gff_file_path}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error parsing GFF file {gff_file_path}: {e}")
        return pd.DataFrame()
        
    df = pd.DataFrame(features)
    
    # Clean up ID column
    if 'ID' not in df.columns and 'locus_tag' in df.columns:
        df['ID'] = df['locus_tag']
    
    return df

def extract_syntenic_region_around_anchor(anchor_gene_name, gff_df_of_target_strain, target_strain_name, neighbors_to_extract=9):
    """Extract syntenic region around an anchor gene with improved ID handling."""
    if gff_df_of_target_strain.empty:
        return None

    # Ensure 'ID', 'locus_tag', 'gene' are string type and stripped
    for col in ['ID', 'locus_tag', 'gene']: 
        if col in gff_df_of_target_strain.columns:
            gff_df_of_target_strain[col] = gff_df_of_target_strain[col].astype(str).str.strip('"').str.strip()

    # Create a unique identifier for each gene based on available columns
    if 'ID' not in gff_df_of_target_strain.columns:
        if 'locus_tag' in gff_df_of_target_strain.columns:
            gff_df_of_target_strain['ID'] = gff_df_of_target_strain['locus_tag']
        else:
            # Generate IDs using strain name, gene name and position
            gff_df_of_target_strain['ID'] = gff_df_of_target_strain.apply(
                lambda row: f"{target_strain_name}_{row.get('gene', 'unknown')}_{row.name}",
                axis=1
            )

    # Convert 'nan' strings to actual NaN values
    gff_df_of_target_strain['ID'] = gff_df_of_target_strain['ID'].replace('nan', pd.NA)
    
    # Fill any remaining NaN IDs with generated IDs
    gff_df_of_target_strain['ID'] = gff_df_of_target_strain['ID'].fillna(
        gff_df_of_target_strain.apply(
            lambda row: f"{target_strain_name}_pos_{row.name}_{row.get('gene', 'unknown')}", 
            axis=1
        )
    )

    # Find anchor gene matches
    anchor_matches = gff_df_of_target_strain[
        gff_df_of_target_strain['gene'].fillna('').str.contains(
            str(anchor_gene_name), regex=False, na=False, case=False
        )
    ]

    anchor_cds_index = -1
    if not anchor_matches.empty:
        for index, row in anchor_matches.iterrows():
            if row.get('type') == 'CDS': # [cite: 79]
                anchor_cds_index = index
                break
    
    if anchor_cds_index == -1: # [cite: 79]
        return None

    # Upstream genes
    upstream_genes_list = []
    current_idx = anchor_cds_index - 1
    while len(upstream_genes_list) < neighbors_to_extract and current_idx >= 0:
        if current_idx < len(gff_df_of_target_strain) and gff_df_of_target_strain.iloc[current_idx].get('type') == 'CDS': #
            upstream_genes_list.append(gff_df_of_target_strain.iloc[current_idx])
        current_idx -= 1
    upstream_df = pd.DataFrame(upstream_genes_list[::-1]) # [cite: 80]

    # Downstream genes
    downstream_genes_list = []
    current_idx = anchor_cds_index + 1
    while len(downstream_genes_list) < neighbors_to_extract and current_idx < len(gff_df_of_target_strain):
        if gff_df_of_target_strain.iloc[current_idx].get('type') == 'CDS': # [cite: 80]
            downstream_genes_list.append(gff_df_of_target_strain.iloc[current_idx])
        current_idx += 1
    downstream_df = pd.DataFrame(downstream_genes_list) # [cite: 80]
    
    # Self (anchor) gene - now uses the definitive 'ID' column
    self_anchor_cols = ['ID', 'type', 'gene', 'start', 'end', 'strand']  # Uses 'ID' directly
    self_anchor_series_data = {}
    for col_loop_var_anchor in self_anchor_cols: 
        self_anchor_series_data[col_loop_var_anchor] = gff_df_of_target_strain.loc[anchor_cds_index, col_loop_var_anchor] if col_loop_var_anchor in gff_df_of_target_strain.columns else None
    self_anchor_df = pd.DataFrame([self_anchor_series_data]) # [cite: 81]
    # No rename of 'Processed_ID' to 'ID' needed here as 'ID' is directly used
    self_anchor_df['position'] = 0 # [cite: 81]

    # Add position; 'ID' column is inherited correctly, no rename needed
    if not upstream_df.empty:
        upstream_df['position'] = range(-len(upstream_df), 0) # [cite: 81]
    if not downstream_df.empty:
        downstream_df['position'] = range(1, len(downstream_df) + 1) # [cite: 82]
    
    extracted_block_df = pd.concat([upstream_df, self_anchor_df, downstream_df], ignore_index=True) # [cite: 82]
    
    output_cols = ['position', 'gene', 'ID', 'type', 'start', 'end', 'strand'] # [cite: 82]
    for col_loop_var_output in output_cols: 
        if col_loop_var_output not in extracted_block_df.columns: extracted_block_df[col_loop_var_output] = pd.NA # [cite: 82]
            
    extracted_block_df['ID'] = extracted_block_df['ID'].fillna('Missing_GFF_ID') # [cite: 82]
    extracted_block_df['gene'] = extracted_block_df['gene'].fillna('Missing_GeneName') # [cite: 82]
    extracted_block_df = extracted_block_df[output_cols] # [cite: 83]
    extracted_block_df['strain'] = target_strain_name # [cite: 83]
    return extracted_block_df

def get_gene_sets_from_context(context_df, gene_col='gene', pos_col='position'): # [cite: 83]
    upstream_set = set()
    downstream_set = set()
    if context_df is None or context_df.empty or gene_col not in context_df.columns or pos_col not in context_df.columns: # [cite: 83]
        return upstream_set, downstream_set
    
    for _, row in context_df.iterrows():
        gene_name = row.get(gene_col) # [cite: 83]
        position_val = row.get(pos_col) # [cite: 83]
        if pd.notna(gene_name) and isinstance(gene_name, str) and pd.notna(position_val): # [cite: 84]
            if position_val < 0: upstream_set.add(gene_name) # [cite: 84]
            elif position_val > 0: downstream_set.add(gene_name) # [cite: 84]
    return upstream_set, downstream_set

def get_anchor_gene_name(reference_context_df, target_position, gene_col='gene', pos_col='position'): # [cite: 84]
    if reference_context_df is None or reference_context_df.empty or \
       pos_col not in reference_context_df.columns or gene_col not in reference_context_df.columns: # [cite: 84]
        return None
    target_row = reference_context_df[reference_context_df[pos_col] == target_position] # [cite: 84]
    if target_row.empty: return None # [cite: 84]
    gene_name = target_row[gene_col].iloc[0] # [cite: 85]
    return gene_name if pd.notna(gene_name) and isinstance(gene_name, str) and gene_name.strip() != '' else None # [cite: 85]

def get_gene_at_position(df, position):
    """Get gene name at a specific position in the DataFrame."""
    if df.loc[df['position'] == position].empty:
        return None
    gene = df.loc[df['position'] == position, 'gene'].fillna('').iloc[0]
    if not isinstance(gene, str) or gene.strip() == '':
        return None
    return gene

def get_upstream_downstream_set(df, column_name='position'):
    """Extract sets of upstream and downstream genes."""
    upstream_genes = set()
    downstream_genes = set()
    
    for _, row in df.iterrows():
        if row[column_name] < 0:
            upstream_genes.add(row['gene'])
        elif row[column_name] > 0:
            downstream_genes.add(row['gene'])
    return upstream_genes, downstream_genes

def check_matches(upstream_genes, downstream_genes, negative_genes, positive_genes):
    """Check matches between gene sets and determine orientation."""
    downstream_matches = 0
    upstream_matches = 0

    neg_downstream_match_count = len(negative_genes.intersection(downstream_genes))
    neg_upstream_match_count = len(negative_genes.intersection(upstream_genes))
    pos_downstream_match_count = len(positive_genes.intersection(downstream_genes))
    pos_upstream_match_count = len(positive_genes.intersection(upstream_genes))

    if neg_downstream_match_count > pos_downstream_match_count:
        downstream_matches = neg_downstream_match_count
        downstream_winner = "negative"
    else:
        downstream_matches = pos_downstream_match_count
        downstream_winner = "positive"
    
    if neg_upstream_match_count > pos_upstream_match_count:
        upstream_matches = neg_upstream_match_count
        upstream_winner = "negative"
    else:
        upstream_matches = pos_upstream_match_count
        upstream_winner = "positive"
        
    return upstream_matches, downstream_matches, upstream_winner, downstream_winner

def extract_data_upstream_downstream(gene, gff_df, strain, neighbors_count=9):
    """Extract upstream and downstream genes around a given gene."""
    if 'ID' in gff_df.columns:
        gff_df['ID'] = gff_df['ID'].str.strip('"').str.strip()
    if 'locus_tag' in gff_df.columns:
        gff_df['locus_tag'] = gff_df['locus_tag'].str.strip('"').str.strip()
    if 'gene' in gff_df.columns:
        gff_df['gene'] = gff_df['gene'].str.strip('"').str.strip()

    if 'locus_tag' in gff_df.columns:
        gff_df['ID'] = gff_df['ID'].fillna(gff_df['locus_tag'])

    gene_matches = gff_df[gff_df['gene'].fillna('').str.contains(str(gene), na=False, case=False)]

    if not gene_matches.empty:
        target_index = gene_matches.index[0]
        for index, row in gene_matches.iterrows():
            if row['type'] == 'CDS':
                target_index = index
                break

        # Extract upstream genes
        upstream_genes = []
        current_idx = target_index - 1
        while len(upstream_genes) < neighbors_count and current_idx >= 0:
            if gff_df.iloc[current_idx]['type'] == 'CDS':
                upstream_genes.append(gff_df.iloc[current_idx])
            current_idx -= 1
        upstream_df = pd.DataFrame(upstream_genes[::-1])
        if not upstream_df.empty:
            upstream_df['position'] = range(-len(upstream_df), 0)

        # Extract downstream genes
        downstream_genes = []
        current_idx = target_index + 1
        while len(downstream_genes) < neighbors_count and current_idx < len(gff_df):
            if gff_df.iloc[current_idx]['type'] == 'CDS':
                downstream_genes.append(gff_df.iloc[current_idx])
            current_idx += 1
        downstream_df = pd.DataFrame(downstream_genes)
        if not downstream_df.empty:
            downstream_df['position'] = range(1, len(downstream_genes) + 1)

        # Process self gene
        self_gene = gff_df.iloc[target_index][['ID', 'type', 'gene', 'start', 'end', 'strand']].copy()
        self_gene['position'] = 0
        self_gene_df = pd.DataFrame([self_gene])

        # Combine all genes
        combined_genes = pd.concat([upstream_df, self_gene_df, downstream_df], ignore_index=True)
        combined_genes['ID'] = combined_genes['ID'].fillna('Missing_ID')
        combined_genes = combined_genes[['position', 'gene', 'ID', 'type', 'start', 'end', 'strand']]
        combined_genes['strain'] = strain

        return combined_genes
    return None

def main():
    parser = argparse.ArgumentParser(description="Finds syntenic regions for soft-core genes reported as missing in some strains.") # [cite: 85]
    parser.add_argument('--input_context_csv', required=True, help="Input data_upstream_downstream.csv from script 1.") # [cite: 85]
    parser.add_argument('--gff_dir', required=True, help="Directory containing GFF files for all strains.") # [cite: 85]
    parser.add_argument('--output_compare_csv', default='compare_genes.csv', help="Output CSV for syntenic regions found (default: compare_genes.csv).") # [cite: 85]
    parser.add_argument('--match_threshold', type=int, default=2, help="Minimum number of matching genes for synteny (default: 2).") #
    args = parser.parse_args()

    try:
        df_contexts = pd.read_csv(args.input_context_csv) # [cite: 86]
        if df_contexts.empty: # [cite: 86]
            print(f"Warning: Input context file '{args.input_context_csv}' is empty. Outputs will be empty.") #
            pd.DataFrame(columns=['position', 'gene', 'ID', 'type', 'start', 'end', 'strand', 'strain', 'central_gene']).to_csv(args.output_compare_csv, index=False) # [cite: 87]
            return
    except FileNotFoundError:
        print(f"Error: Input context file '{args.input_context_csv}' not found.") # [cite: 87]
        return
    except Exception as e:
        print(f"Error reading {args.input_context_csv}: {e}") #
        return

    hardcoded_strain_list = ["CFT073","ECOR11_PROKKA","ECOR14_PROKKA","ECOR16_PROKKA","ECOR23_PROKKA","ECOR44_PROKKA","ECOR47_PROKKA","ECOR48_PROKKA","ECOR49_PROKKA","ECOR4_PROKKA","ECOR50_PROKKA","ECOR51_PROKKA","ECOR52_PROKKA","ECOR53_PROKKA","ECOR54_PROKKA","ECOR55_PROKKA","ECOR56_PROKKA","ECOR60_PROKKA","ECOR61_PROKKA","ECOR62_PROKKA","ECOR64_PROKKA","ECOR66_PROKKA","ECOR6_PROKKA","ECOR8_PROKKA","KE21","KE24","KE25","KE26","KE4","KE40","KE41","KE46","KE47","KE48_PROKKA","KE50","KE51","KE54","KE55","KE58","LRPF007","LRPF62","MG1655","UTI89","W3110"] # [cite: 88]
    all_possible_strains = set(hardcoded_strain_list) # [cite: 88]
    if 'strain' in df_contexts.columns: # [cite: 88]
        all_possible_strains.update(df_contexts['strain'].unique()) # [cite: 88]

    missing_strains_map = {} # [cite: 88]
    if 'central_gene' in df_contexts.columns and 'strain' in df_contexts.columns: # [cite: 88]
        for central_gene_name, group_data in df_contexts.groupby('central_gene'): # [cite: 88]
            present_strains_for_gene = set(group_data['strain'].unique()) # [cite: 88]
            missing_strains_map[central_gene_name] = list(all_possible_strains - present_strains_for_gene) # [cite: 88]
    else:
        print("Error: 'central_gene' or 'strain' columns not found in input context CSV. Cannot proceed.") #
        return

    gff_data_cache = {} # [cite: 90]
    print("Pre-loading GFF files for synteny matching...") # [cite: 90]
    strains_needing_gff = set(ms for cg_miss_list in missing_strains_map.values() for ms in cg_miss_list if isinstance(cg_miss_list, list)) # [cite: 90]
    for strain_name_iter in strains_needing_gff: # [cite: 90]
        gff_file_path = os.path.join(args.gff_dir, f"{strain_name_iter}.gff") # [cite: 90]
        if os.path.exists(gff_file_path) and strain_name_iter not in gff_data_cache: # [cite: 90]
            parsed_gff = parse_gff_to_df(gff_file_path) # [cite: 90]
            if not parsed_gff.empty: # [cite: 91]
                gff_data_cache[strain_name_iter] = parsed_gff # [cite: 91]
    print(f"Finished pre-loading {len(gff_data_cache)} GFF files.") # [cite: 91]

    output_compare_cols = ['position', 'gene', 'ID', 'type', 'start', 'end', 'strand', 'strain', 'central_gene'] # [cite: 91]

    synteny_found_records_dfs = [] # [cite: 91]
    synteny_block_written_for_pair = set()  # [cite: 91]

    if 'central_gene' in df_contexts.columns and 'strain' in df_contexts.columns: # [cite: 91]
        for (ref_roary_gene, ref_strain), reference_context_df in df_contexts.groupby(['central_gene', 'strain']): # [cite: 91]
            ref_upstream_set, ref_downstream_set = get_gene_sets_from_context(reference_context_df) # [cite: 92]
            strains_where_ref_roary_gene_is_missing = missing_strains_map.get(ref_roary_gene, []) # [cite: 92]

            for missing_strain_name in strains_where_ref_roary_gene_is_missing: # [cite: 92]
                if missing_strain_name not in gff_data_cache: # [cite: 92]
                    continue
                gff_df_for_missing_strain = gff_data_cache[missing_strain_name] # [cite: 93
                found_match_for_this_missing_strain = False # [cite: 93]

                for anchor_search_direction in ["upstream", "downstream"]: # [cite: 93]
                    if found_match_for_this_missing_strain: break # [cite: 93]
                    ref_pos_range = range(-5, 0) if anchor_search_direction == "upstream" else range(1, 6) # [cite: 93]
                    
                    for ref_pos in ref_pos_range: # [cite: 94]
                        anchor_gene_name_in_ref = get_anchor_gene_name(reference_context_df, ref_pos) # [cite: 94]
                        if not anchor_gene_name_in_ref: # [cite: 94]
                            continue
                        
                        # Pass a copy of the gff_df to avoid modifying the cache
                        extracted_block_in_missing_strain = extract_syntenic_region_around_anchor(
                            anchor_gene_name_in_ref, gff_df_for_missing_strain.copy(), missing_strain_name, neighbors_to_extract=9 
                        ) #

                        if extracted_block_in_missing_strain is not None and not extracted_block_in_missing_strain.empty: # [cite: 96]
                            block_upstream_set, block_downstream_set = get_gene_sets_from_context(extracted_block_in_missing_strain) # [cite: 96]
                            num_upstream_matches = len(ref_upstream_set.intersection(block_upstream_set)) # [cite: 97
                            num_downstream_matches = len(ref_downstream_set.intersection(block_downstream_set)) # [cite: 97

                            if num_upstream_matches >= args.match_threshold and num_downstream_matches >= args.match_threshold: # [cite: 97]
                                if (ref_roary_gene, missing_strain_name) not in synteny_block_written_for_pair: # [cite: 100]
                                    extracted_block_in_missing_strain['central_gene'] = ref_roary_gene # [cite: 101]
                                    # Ensure IDs are present in the extracted block
                                    if 'ID' not in extracted_block_in_missing_strain.columns:
                                        extracted_block_in_missing_strain['ID'] = extracted_block_in_missing_strain.index.map(lambda x: f"{missing_strain_name}_gene_{x}")
                                    synteny_found_records_dfs.append(extracted_block_in_missing_strain) # [cite: 101]
                                    synteny_block_written_for_pair.add((ref_roary_gene, missing_strain_name)) # [cite: 101]
                                found_match_for_this_missing_strain = True # [cite: 102]
                                break # [cite: 102]
    
    compare_genes_output_df = pd.DataFrame(columns=output_compare_cols) # [cite: 102]
    if synteny_found_records_dfs: # [cite: 102]
        compare_genes_output_df = pd.concat(synteny_found_records_dfs, ignore_index=True) # [cite: 102]
        # Ensure all expected columns are present, fill with NA if not, and set order
        for col_loop_var_final_df in output_compare_cols: 
            if col_loop_var_final_df not in compare_genes_output_df.columns:
                compare_genes_output_df[col_loop_var_final_df] = pd.NA
        compare_genes_output_df = compare_genes_output_df[output_compare_cols] # [cite: 102]
        
    # Add directory creation before writing output
    output_dir = os.path.dirname(args.output_compare_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    try:
        compare_genes_output_df.to_csv(args.output_compare_csv, index=False)
        print(f"Syntenic regions data ({len(compare_genes_output_df)} rows) written to {args.output_compare_csv}")
    except PermissionError:
        alt_output = "compare_genes.csv"  # Fallback to current directory
        compare_genes_output_df.to_csv(alt_output, index=False)
        print(f"Permission denied for {args.output_compare_csv}")
        print(f"Writing output to current directory instead: {alt_output}")
    except Exception as e:
        print(f"Error writing output: {e}")
        return

if __name__ == '__main__':
    main()