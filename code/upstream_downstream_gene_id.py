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
    if gff_df.empty:
        return None

    for col_name in ['ID', 'locus_tag', 'gene']:
        if col_name in gff_df.columns:
            gff_df[col_name] = gff_df[col_name].astype(str).str.strip('"').str.strip()

    if 'ID' in gff_df.columns and 'locus_tag' in gff_df.columns:
        gff_df['Processed_ID'] = gff_df['ID'].fillna(gff_df['locus_tag'])
    elif 'ID' in gff_df.columns:
        gff_df['Processed_ID'] = gff_df['ID']
    elif 'locus_tag' in gff_df.columns:
        gff_df['Processed_ID'] = gff_df['locus_tag']
    else:
        return None
    
    id_matches = gff_df[gff_df['Processed_ID'].fillna('') == str(target_gene_id_from_roary)]
    
    target_cds_index = -1
    if not id_matches.empty:
        for index, row in id_matches.iterrows():
            if row.get('type') == 'CDS':
                target_cds_index = index
                break
    
    if target_cds_index == -1:
        return None

    upstream_genes_list = []
    current_idx = target_cds_index - 1
    while len(upstream_genes_list) < neighbors_count and current_idx >= 0:
        if current_idx < len(gff_df) and gff_df.iloc[current_idx].get('type') == 'CDS':
            upstream_genes_list.append(gff_df.iloc[current_idx])
        current_idx -= 1
    upstream_df = pd.DataFrame(upstream_genes_list).iloc[::-1]

    downstream_genes_list = []
    current_idx = target_cds_index + 1
    while len(downstream_genes_list) < neighbors_count and current_idx < len(gff_df):
        if gff_df.iloc[current_idx].get('type') == 'CDS':
            downstream_genes_list.append(gff_df.iloc[current_idx])
        current_idx += 1
    downstream_df = pd.DataFrame(downstream_genes_list)

    self_gene_cols_to_extract = ['ID', 'type', 'gene', 'start', 'end', 'strand']
    self_gene_series_data = {}
    for col in self_gene_cols_to_extract:
        self_gene_series_data[col] = gff_df.loc[target_cds_index, col] if col in gff_df.columns else None
    self_gene_series_data['ID'] = target_gene_id_from_roary

    self_gene_df = pd.DataFrame([self_gene_series_data])
    self_gene_df['Position'] = 0

    if not upstream_df.empty:
        upstream_df['Position'] = range(-len(upstream_df), 0)
        if 'Processed_ID' in upstream_df.columns:
            upstream_df.rename(columns={'Processed_ID': 'ID'}, inplace=True)
        elif 'ID' not in upstream_df.columns and 'locus_tag' in upstream_df.columns:
             upstream_df.rename(columns={'locus_tag': 'ID'}, inplace=True)

    if not downstream_df.empty:
        downstream_df['Position'] = range(1, len(downstream_df) + 1)
        if 'Processed_ID' in downstream_df.columns:
            downstream_df.rename(columns={'Processed_ID': 'ID'}, inplace=True)
        elif 'ID' not in downstream_df.columns and 'locus_tag' in downstream_df.columns:
             downstream_df.rename(columns={'locus_tag': 'ID'}, inplace=True)
            
    combined_df = pd.concat([upstream_df, self_gene_df, downstream_df], ignore_index=True)
    
    output_cols = ['Position', 'gene', 'ID', 'type', 'start', 'end', 'strand']
    for col in output_cols:
        if col not in combined_df.columns:
            combined_df[col] = pd.NA
    
    combined_df['ID'] = combined_df['ID'].fillna('Missing_GFF_ID_In_Neighbor')
    combined_df['gene'] = combined_df['gene'].fillna('Missing_GeneName')
    
    combined_df = combined_df[output_cols]
    combined_df['strain'] = strain_name
    return combined_df

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