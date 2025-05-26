# upstream_downstream_gene_id.py
import os
import pandas as pd
import csv
import argparse
import re

def parse_gff_to_df(gff_file_path):
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
                            attr_dict[key] = value.strip()
                    feature = {
                        'seqid': seqid, 'source': source, 'type': ftype,
                        'start': int(start), 'end': int(end), 'score': score,
                        'strand': strand, 'phase': phase, **attr_dict
                    }
                    features.append(feature)
    except FileNotFoundError:
        print(f"Warning: GFF file not found at {gff_file_path}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Warning: Error parsing GFF file {gff_file_path}: {e}")
        return pd.DataFrame()
    return pd.DataFrame(features)

def extract_data_upstream_downstream(target_gene_id_from_roary, gff_df, strain_name, neighbors_count=5):
    """Extract upstream and downstream genes around a given Roary ID."""
    # Reset index to ensure unique indices
    gff_df = gff_df.reset_index(drop=True)
    
    # Sanitize string columns
    if 'ID' in gff_df.columns:
        gff_df['ID'] = gff_df['ID'].str.strip('"').str.strip()
    if 'locus_tag' in gff_df.columns:
        gff_df['locus_tag'] = gff_df['locus_tag'].str.strip('"').str.strip()
    if 'gene' in gff_df.columns:
        gff_df['gene'] = gff_df['gene'].str.strip('"').str.strip()

    # Merge ID and locus_tag into a single ID column
    if 'locus_tag' in gff_df.columns:
        gff_df['ID'] = gff_df['ID'].fillna(gff_df['locus_tag'])

    # Use loose matching to find the ID in the ID column
    id_matches = gff_df[gff_df['ID'].fillna('').str.contains(str(target_gene_id_from_roary), na=False, case=False)]

    if not id_matches.empty:
        target_index = id_matches.index[0]  # Get the index of the target gene
        for index, row in id_matches.iterrows():
            if row['type'] == 'CDS':
                target_index = index
                break

        # --- Upstream genes ---
        upstream_genes = []
        upstream_genes_index = target_index - 1
        while len(upstream_genes) < neighbors_count and upstream_genes_index >= 0:
            if gff_df.iloc[upstream_genes_index]['type'] == 'CDS':
                upstream_genes.append(gff_df.iloc[upstream_genes_index])
            upstream_genes_index -= 1

        upstream_genes = pd.DataFrame(upstream_genes).iloc[::-1]

        # --- Downstream genes ---
        downstream_genes = []
        downstream_genes_index = target_index + 1
        while len(downstream_genes) < neighbors_count and downstream_genes_index < len(gff_df):
            if gff_df.iloc[downstream_genes_index]['type'] == 'CDS':
                downstream_genes.append(gff_df.iloc[downstream_genes_index])
            downstream_genes_index += 1

        downstream_genes = pd.DataFrame(downstream_genes)

        # Assign positions (update case to match output header)
        if not upstream_genes.empty:
            upstream_genes['position'] = range(-len(upstream_genes), 0)
        if not downstream_genes.empty:
            downstream_genes['position'] = range(1, len(downstream_genes) + 1)

        # Extract and prepare self gene
        self_gene = gff_df.iloc[target_index][['ID', 'type', 'gene', 'start', 'end', 'strand']].copy()
        self_gene['position'] = 0  # Changed Position to position
        self_gene_df = pd.DataFrame([self_gene])

        # Combine all genes
        combined_genes = pd.concat([upstream_genes, self_gene_df, downstream_genes], ignore_index=True)
        
        # Fill missing values and select columns (note lowercase 'position')
        combined_genes['ID'] = combined_genes['ID'].fillna('Missing_ID')
        combined_genes = combined_genes[['position', 'gene', 'ID', 'type', 'start', 'end', 'strand']]
        combined_genes['strain'] = strain_name

        return combined_genes
    else:
        print(f"No row found with ID: {target_gene_id_from_roary} for strain: {strain_name}")
        return None
def add_core_status_and_central_gene(combined_genes_df, roary_gene_info_map, central_roary_gene_name_for_this_block):
    if combined_genes_df is None or combined_genes_df.empty:
        return None
    
    combined_genes_df['core_status'] = 'unknown' 
    combined_genes_df['isolate_count'] = pd.NA
    combined_genes_df['central_gene'] = central_roary_gene_name_for_this_block

    for index, row in combined_genes_df.iterrows():
        gene_id_for_lookup = row.get('ID') 
        
        if gene_id_for_lookup in roary_gene_info_map:
            roary_name, isolate_count, core_status = roary_gene_info_map[gene_id_for_lookup]
            combined_genes_df.loc[index, 'core_status'] = core_status
            combined_genes_df.loc[index, 'isolate_count'] = isolate_count
            current_gene_name = row.get('gene')
            if pd.isna(current_gene_name) or current_gene_name == 'Missing_GeneName':
                 combined_genes_df.loc[index, 'gene'] = roary_name
    return combined_genes_df

def main():
    parser = argparse.ArgumentParser(description="Extracts upstream/downstream gene context for specific soft-core genes.")
    parser.add_argument('--presence_absence_csv', required=True, help="Input Roary gene_presence_absence.refmt.csv file.")
    parser.add_argument('--gff_dir', required=True, help="Directory containing GFF files for all strains.")
    parser.add_argument('--missing_gene_out', default='missing_gene_data.csv', help="Output CSV file for missing gene data (default: missing_gene_data.csv).")
    parser.add_argument('--upstream_downstream_out', default='data_upstream_downstream.csv', help="Output CSV for upstream/downstream gene data (default: data_upstream_downstream.csv).")
    args = parser.parse_args()

    try:
        df_roary = pd.read_csv(args.presence_absence_csv)
    except FileNotFoundError:
        print(f"Error: Input file '{args.presence_absence_csv}' not found.")
        return
    except Exception as e:
        print(f"Error reading {args.presence_absence_csv}: {e}")
        return

    if 'Gene' in df_roary.columns and 'Non-unique Gene name' in df_roary.columns:
        df_roary['Gene'] = df_roary.apply(
            lambda r: r['Non-unique Gene name'] if pd.notna(r['Non-unique Gene name']) and r['Non-unique Gene name'].strip() != '' else r['Gene'],
            axis=1
        )

    strain_columns = ["CFT073","ECOR11_PROKKA","ECOR14_PROKKA","ECOR16_PROKKA","ECOR23_PROKKA","ECOR44_PROKKA","ECOR47_PROKKA","ECOR48_PROKKA","ECOR49_PROKKA","ECOR4_PROKKA","ECOR50_PROKKA","ECOR51_PROKKA","ECOR52_PROKKA","ECOR53_PROKKA","ECOR54_PROKKA","ECOR55_PROKKA","ECOR56_PROKKA","ECOR60_PROKKA","ECOR61_PROKKA","ECOR62_PROKKA","ECOR64_PROKKA","ECOR66_PROKKA","ECOR6_PROKKA","ECOR8_PROKKA","KE21","KE24","KE25","KE26","KE4","KE40","KE41","KE46","KE47","KE48_PROKKA","KE50","KE51","KE54","KE55","KE58","LRPF007","LRPF62","MG1655","UTI89","W3110"]
    
    roary_gene_info_map = {} 
    for index, row in df_roary.iterrows():
        num_isolates = row.get('No. isolates', 0)
        core_stat = "unknown"
        if num_isolates == 44: core_stat = "core"
        elif num_isolates == 42 or num_isolates == 43: core_stat = "soft_core" 
        elif num_isolates < 42: core_stat = "non_core" 

        roary_cluster_name = row.get('Gene', f"Cluster_{index}")

        for strain_col_name in strain_columns:
            if strain_col_name in row and pd.notna(row[strain_col_name]):
                gene_ids_in_cell = str(row[strain_col_name]).split(';') 
                for single_roary_id in gene_ids_in_cell:
                    single_roary_id = single_roary_id.strip()
                    if single_roary_id and single_roary_id not in roary_gene_info_map:
                        roary_gene_info_map[single_roary_id] = (roary_cluster_name, num_isolates, core_stat)

    df_soft_core_candidates = df_roary[df_roary.get('No. isolates', 0) == 43].copy()

    missing_gene_output_list = []
    if not df_soft_core_candidates.empty:
        for index, row in df_soft_core_candidates.iterrows():
            roary_cluster_name = row.get('Gene', f"Cluster_{index}")
            for strain_col_name in strain_columns:
                if strain_col_name in row and pd.isnull(row[strain_col_name]): 
                    missing_gene_output_list.append({'Gene': roary_cluster_name, 'Strain': strain_col_name})
    
    pd.DataFrame(missing_gene_output_list).to_csv(args.missing_gene_out, index=False)
    print(f"Missing gene data for '43-isolate genes' written to {args.missing_gene_out}")

    gene_contexts_to_extract = []
    if not df_soft_core_candidates.empty:
        for index, row in df_soft_core_candidates.iterrows():
            roary_cluster_name = row.get('Gene', f"Cluster_{index}")
            for strain_col_name in strain_columns:
                if strain_col_name in row and pd.notna(row[strain_col_name]):
                    gene_ids_in_cell = str(row[strain_col_name]).split(';')
                    for single_roary_id in gene_ids_in_cell:
                        single_roary_id = single_roary_id.strip()
                        if single_roary_id:
                            gene_contexts_to_extract.append({
                                'roary_cluster_name': roary_cluster_name,
                                'strain_present_in': strain_col_name,
                                'roary_id_in_strain': single_roary_id
                            })
    
    all_gff_data = {}
    print("Pre-loading GFF files for context extraction...")
    strains_for_gff_loading = set(item['strain_present_in'] for item in gene_contexts_to_extract)
    for s_name in strains_for_gff_loading:
        gff_path = os.path.join(args.gff_dir, f"{s_name}.gff")
        if s_name not in all_gff_data: 
             parsed_gff = parse_gff_to_df(gff_path)
             if not parsed_gff.empty:
                 all_gff_data[s_name] = parsed_gff
    print(f"Finished pre-loading {len(all_gff_data)} GFF files.")

    all_processed_context_dfs = []
    output_header = ['position', 'gene', 'ID', 'type', 'start', 'end', 'strand', 'strain', 'core_status', 'isolate_count', 'central_gene']

    for item in gene_contexts_to_extract:
        strain_name = item['strain_present_in']
        roary_id = item['roary_id_in_strain']
        roary_cluster_name = item['roary_cluster_name']

        if strain_name not in all_gff_data:
            continue
        
        gff_df_for_strain = all_gff_data[strain_name]
        try:
            context_df = extract_data_upstream_downstream(roary_id, gff_df_for_strain.copy(), strain_name)
            if context_df is not None and not context_df.empty:
                context_df_annotated = add_core_status_and_central_gene(context_df, roary_gene_info_map, roary_cluster_name)
                if context_df_annotated is not None and not context_df_annotated.empty:
                    all_processed_context_dfs.append(context_df_annotated)
        except Exception as e_context:
            print(f"Error extracting context for Roary ID {roary_id} in strain {strain_name} (Cluster: {roary_cluster_name}): {e_context}")

    if all_processed_context_dfs:
        final_df = pd.concat(all_processed_context_dfs, ignore_index=True)
        for col_hdr in output_header: 
            if col_hdr not in final_df.columns:
                final_df[col_hdr] = pd.NA 
        final_df = final_df[output_header] 
        final_df.to_csv(args.upstream_downstream_out, index=False, header=True)
        print(f"Upstream/downstream context data ({len(final_df)} rows) written to {args.upstream_downstream_out}")
    else:
        pd.DataFrame(columns=output_header).to_csv(args.upstream_downstream_out, index=False, header=True)
        print(f"No upstream/downstream context data generated. Empty file created: {args.upstream_downstream_out}")

if __name__ == '__main__':
    main()