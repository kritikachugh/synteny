#!/usr/bin/env python3
import pandas as pd
import os
import csv
import argparse

# This list is considered static based on the dataset's structure.
# If it needs to be dynamic, it would have to be derived from the input CSV.
strain_columns = ["CFT073","ECOR11_PROKKA","ECOR14_PROKKA","ECOR16_PROKKA","ECOR23_PROKKA","ECOR44_PROKKA","ECOR47_PROKKA","ECOR48_PROKKA","ECOR49_PROKKA","ECOR4_PROKKA","ECOR50_PROKKA","ECOR51_PROKKA","ECOR52_PROKKA","ECOR53_PROKKA","ECOR54_PROKKA","ECOR55_PROKKA","ECOR56_PROKKA","ECOR60_PROKKA","ECOR61_PROKKA","ECOR62_PROKKA","ECOR64_PROKKA","ECOR66_PROKKA","ECOR6_PROKKA","ECOR8_PROKKA","KE21","KE24","KE25","KE26","KE4","KE40","KE41","KE46","KE47","KE48_PROKKA","KE50","KE51","KE54","KE55","KE58","LRPF007","LRPF62","MG1655","UTI89","W3110"]

def parse_gff_to_df(gff_file):
  """
  Parses a GFF file and returns a pandas DataFrame of features.
  """
  features = []
  with open(gff_file, 'r') as f:
    for line in f:
      if line.startswith('#'):
        continue
      fields = line.strip().split('\t')
      if len(fields) >= 9:
        seqid, source, ftype, start, end, score, strand, phase, attributes = fields[:9]
        attr_dict = {}
        for attr in attributes.split(';'):
          if '=' in attr:
            key, value = attr.split('=', 1) # Split only on the first '='
            attr_dict[key] = value

        feature = {
            'seqid': seqid,
            'source': source,
            'type': ftype,
            'start': int(start),
            'end': int(end),
            'score': score,
            'strand': strand,
            'phase': phase,
            **attr_dict
        }
        features.append(feature)
  return pd.DataFrame(features)


def extract_data_upstream_downstream(gene, gff_df, strain, neighbors_count=5):
    # Sanitize columns to remove extra spaces and quotes
    for col in ['ID', 'locus_tag', 'gene']:
        if col in gff_df.columns:
            gff_df[col] = gff_df[col].astype(str).str.strip('"').str.strip()

    # Merge 'locus_tag' into 'ID'
    if 'locus_tag' in gff_df.columns:
        gff_df['ID'] = gff_df['ID'].fillna(gff_df['locus_tag'])

    # Find the index of the target gene
    gene_matches = gff_df[gff_df['gene'].fillna('').str.fullmatch(gene, case=False)]

    if gene_matches.empty:
        # Fallback to partial match if full match fails
        gene_matches = gff_df[gff_df['gene'].fillna('').str.contains(gene, na=False, case=False)]

    if not gene_matches.empty:
        target_index = -1
        # Prioritize CDS type features
        for index, row in gene_matches.iterrows():
            if row['type'] == 'CDS':
                target_index = index
                break
        if target_index == -1:
            target_index = gene_matches.index[0]

        # --- Upstream genes ---
        upstream_genes = []
        upstream_genes_index = target_index - 1
        while len(upstream_genes) < neighbors_count and upstream_genes_index >= 0:
            if gff_df.iloc[upstream_genes_index]['type'] == 'CDS':
                upstream_genes.append(gff_df.iloc[upstream_genes_index])
            upstream_genes_index -= 1
        
        # Positions are relative to the *found* gene, so from -1 to -neighbors_count
        if upstream_genes:
            upstream_df = pd.DataFrame(upstream_genes[::-1]) # Farthest is first
            upstream_df['position'] = range(-len(upstream_df), 0)
        else:
            upstream_df = pd.DataFrame()


        # --- Downstream genes ---
        downstream_genes = []
        downstream_genes_index = target_index + 1
        while len(downstream_genes) < neighbors_count and downstream_genes_index < len(gff_df):
            if gff_df.iloc[downstream_genes_index]['type'] == 'CDS':
                downstream_genes.append(gff_df.iloc[downstream_genes_index])
            downstream_genes_index += 1
        
        if downstream_genes:
            downstream_df = pd.DataFrame(downstream_genes)
            # Positions are relative to the *found* gene, from 2 to neighbors_count+1
            downstream_df['position'] = range(2, len(downstream_df) + 2)
        else:
            downstream_df = pd.DataFrame()


        # Extract the target gene (self gene)
        self_gene = pd.DataFrame([gff_df.loc[target_index]])
        self_gene['position'] = 1 # Position of the found gene is always 1

        combined_genes = pd.concat([upstream_df, self_gene, downstream_df], ignore_index=True)
        combined_genes['ID'] = combined_genes['ID'].fillna('Missing_ID')
        
        # Ensure all necessary columns exist, fill with None if not
        for col in ['gene', 'type', 'start', 'end', 'strand']:
            if col not in combined_genes.columns:
                combined_genes[col] = None

        combined_genes = combined_genes[['position', 'gene', 'ID', 'type', 'start', 'end', 'strand']]
        combined_genes['strain'] = strain

        return combined_genes
    else:
        # To avoid clutter, we can comment this out or use a logging library
        # print(f"Warning: Gene '{gene}' not found for strain '{strain}'.")
        return None

def get_upstream_downstream_set(df, column_name):
  upstream_genes = set(df[df[column_name] < 0]['gene'].dropna())
  downstream_genes = set(df[df[column_name] > 0]['gene'].dropna())
  return upstream_genes, downstream_genes

def check_matches(upstream_genes, downstream_genes, negative_genes, positive_genes):
  neg_downstream_match_count = len(negative_genes.intersection(downstream_genes))
  pos_downstream_match_count = len(positive_genes.intersection(downstream_genes))
  
  downstream_matches = max(neg_downstream_match_count, pos_downstream_match_count)
  
  neg_upstream_match_count = len(negative_genes.intersection(upstream_genes))
  pos_upstream_match_count = len(positive_genes.intersection(upstream_genes))

  upstream_matches = max(neg_upstream_match_count, pos_upstream_match_count)

  return upstream_matches, downstream_matches

def get_gene_at_position(df, position):
  gene_series = df.loc[df['position'] == position, 'gene']
  if gene_series.empty:
    return None
  gene = gene_series.fillna('').iloc[0]
  if not isinstance(gene, str) or gene.strip() == '':
    return None
  return gene

def get_full_gff_file_path(gff_dir, strain):
  return os.path.join(gff_dir, strain + ".gff")


def main(args):
    """Main execution function"""
    # Read input files and configurations from args
    input_csv_path = args.input_context_csv
    gff_dir = args.gff_dir
    output_csv_path = args.output_compare_csv
    threshold_match = args.match_threshold

    # Derive the matches file name from the main output file name
    base, ext = os.path.splitext(output_csv_path)
    matches_csv_path = f"{base}_matches{ext}"
    
    print(f"Reading input context from: {input_csv_path}")
    upstream_downstream_df = pd.read_csv(input_csv_path)

    # Group by the central gene to find missing strains for each
    missing_strains_dict = {}
    for central_gene, gene_data in upstream_downstream_df.groupby('central_gene'):
        present_strains = set(gene_data['strain'])
        missing_strains = set(strain_columns) - present_strains
        missing_strains_dict[central_gene] = list(missing_strains)

    print(f"Loading GFF files from directory: {gff_dir}")
    gff_dfs = {}
    for strain_column in strain_columns:
        gff_file_path = get_full_gff_file_path(gff_dir, strain_column)
        if os.path.exists(gff_file_path):
            gff_dfs[strain_column] = parse_gff_to_df(gff_file_path)
        else:
            print(f"Warning: GFF file not found for strain {strain_column} at {gff_file_path}")

    # Prepare output files
    csv_columns = ['position', 'gene', 'ID', 'type', 'start', 'end', 'strand', 'strain', 'central_gene']
    matches_csv_columns = ['central_gene', 'missing_strain', 'reference_strain', 'downstream_matches', 'upstream_matches']
    matches_df_list = []

    print(f"Starting comparison... Writing synteny results to {output_csv_path}")
    seen_matches = set()
    with open(output_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(csv_columns)

        # Iterate through each known context (central gene in a specific strain)
        for (central_gene, strain), group in upstream_downstream_df.groupby(['central_gene', 'strain']):
            filtered_df = group.copy()
            
            # Get the upstream/downstream gene sets for the known context
            upstream_genes, downstream_genes = get_upstream_downstream_set(filtered_df, 'position')
            
            # Get the list of strains where this central gene is missing
            missing_strains = missing_strains_dict.get(central_gene, [])

            for missing_strain in missing_strains:
                if missing_strain not in gff_dfs:
                    continue # Skip if GFF for missing strain was not loaded

                gff_df = gff_dfs[missing_strain]
                found_match = False
                
                # Search using genes upstream of the known context
                for position in range(-1, -6, -1):
                    anchor_gene = get_gene_at_position(filtered_df, position)
                    if anchor_gene is None:
                        continue
                    
                    result = extract_data_upstream_downstream(anchor_gene, gff_df.copy(), missing_strain, neighbors_count=9)
                    
                    if result is not None:
                        found_match = True
                        result['central_gene'] = central_gene
                        
                        # Compare the newly found context with the original known context
                        new_up_genes, new_down_genes = get_upstream_downstream_set(result, 'position')
                        upstream_matches, downstream_matches = check_matches(upstream_genes, downstream_genes, new_up_genes, new_down_genes)
                        
                        if upstream_matches >= threshold_match and downstream_matches >= threshold_match:
                            match_key = (central_gene, missing_strain)
                            if match_key not in seen_matches:
                                seen_matches.add(match_key)
                                result.to_csv(csvfile, header=False, index=False, columns=csv_columns)
                            
                            matches_df_list.append([central_gene, missing_strain, strain, downstream_matches, upstream_matches])
                        break # Move to the next missing strain once a match is processed
                
                if found_match:
                    continue

                # If no match found upstream, search using genes downstream
                for position in range(2, 6): # Original central gene is at pos 1, so start from 2
                    anchor_gene = get_gene_at_position(filtered_df, position)
                    if anchor_gene is None:
                        continue

                    result = extract_data_upstream_downstream(anchor_gene, gff_df.copy(), missing_strain, neighbors_count=9)
                    
                    if result is not None:
                        result['central_gene'] = central_gene
                        new_up_genes, new_down_genes = get_upstream_downstream_set(result, 'position')
                        upstream_matches, downstream_matches = check_matches(upstream_genes, downstream_genes, new_up_genes, new_down_genes)

                        if upstream_matches >= threshold_match and downstream_matches >= threshold_match:
                            match_key = (central_gene, missing_strain)
                            if match_key not in seen_matches:
                                seen_matches.add(match_key)
                                result.to_csv(csvfile, header=False, index=False, columns=csv_columns)

                            matches_df_list.append([central_gene, missing_strain, strain, downstream_matches, upstream_matches])
                        break # Move to the next missing strain
    
    # Save the summary of matches
    if matches_df_list:
        matches_df = pd.DataFrame(matches_df_list, columns=matches_csv_columns)
        matches_df.to_csv(matches_csv_path, index=False)
        print(f"Wrote match summary to {matches_csv_path}")
    
    print("Processing complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find and compare syntenic regions for genes missing in certain strains. "
                    "It uses a known gene context from one strain to find a similar region in another "
                    "strain where the central gene is absent."
    )
    parser.add_argument(
        "--input_context_csv",
        type=str,
        required=True,
        help="Path to the input CSV file containing known gene contexts (e.g., data_upstream_downstream.csv)."
    )
    parser.add_argument(
        "--gff_dir",
        type=str,
        required=True,
        help="Path to the directory containing all GFF files, named as <strain_name>.gff."
    )
    parser.add_argument(
        "--output_compare_csv",
        type=str,
        required=True,
        help="Path for the output CSV file that will contain the detailed syntenic regions found."
    )
    parser.add_argument(
        "--match_threshold",
        type=int,
        default=2,
        help="The minimum number of matching genes required in both upstream and downstream regions to be considered a syntenic match. Default is 2."
    )
    
    args = parser.parse_args()
    main(args)