# cds_sequence.py
import os
import pandas as pd
import csv
import argparse

def get_cds_file_path(strain_name_val, cds_files_directory):
  full_path = os.path.join(cds_files_directory, strain_name_val + ".CDS.txt")
  return full_path

def get_sequence_by_gene_id_from_df(target_gene_id, cds_df_parsed):
    if cds_df_parsed is None or cds_df_parsed.empty or 'gene_id' not in cds_df_parsed or 'sequence' not in cds_df_parsed:
        return None
    try:
        gene_data_row = cds_df_parsed[cds_df_parsed['gene_id'].astype(str) == str(target_gene_id)]
        if not gene_data_row.empty:
            sequence_val = gene_data_row['sequence'].iloc[0]
            return sequence_val
        return None
    except Exception:
        return None

def read_cds_file_to_df(file_path_val):
  try:
    df_cds = pd.read_csv(file_path_val, sep='\t', header=None, dtype=str, keep_default_na=False, na_values=[''])
    if df_cds.shape[1] < 8:
        return None
    cds_headers_list = ['accession_no', 'start', 'end', 'gene_id', 'na_col', 'orientation', 'gene_name', 'sequence']
    df_cds.columns = cds_headers_list[:df_cds.shape[1]]
    return df_cds
  except FileNotFoundError:
    return None
  except pd.errors.EmptyDataError:
    return None
  except Exception as e:
    print(f"Error reading or parsing CDS file {file_path_val}: {e}")
    return None

def main():
    parser = argparse.ArgumentParser(description="Extract CDS sequences based on gene presence/absence data.")
    parser.add_argument('--presence_absence_csv', required=True, help="Path to gene_presence_absence.refmt.csv")
    parser.add_argument('--cds_dir', required=True, help="Directory containing CDS files (e.g., <strain_name>.CDS.txt).")
    parser.add_argument('--output_csv', required=True, help="Output path for sequence_data.csv.")
    args = parser.parse_args()

    try:
        df_roary = pd.read_csv(args.presence_absence_csv)
    except FileNotFoundError:
        print(f"Error: File not found {args.presence_absence_csv}")
        return
    except Exception as e:
        print(f"Error reading {args.presence_absence_csv}: {e}")
        return

    if 'Gene' in df_roary.columns and 'Non-unique Gene name' in df_roary.columns:
        df_roary['Gene'] = df_roary.apply(
            lambda r: r['Non-unique Gene name'] if pd.notna(r['Non-unique Gene name']) and r['Non-unique Gene name'].strip() != '' else r['Gene'],
            axis=1
        )

    strain_names_list = ["CFT073","ECOR11_PROKKA","ECOR14_PROKKA","ECOR16_PROKKA","ECOR23_PROKKA","ECOR44_PROKKA","ECOR47_PROKKA","ECOR48_PROKKA","ECOR49_PROKKA","ECOR4_PROKKA","ECOR50_PROKKA","ECOR51_PROKKA","ECOR52_PROKKA","ECOR53_PROKKA","ECOR54_PROKKA","ECOR55_PROKKA","ECOR56_PROKKA","ECOR60_PROKKA","ECOR61_PROKKA","ECOR62_PROKKA","ECOR64_PROKKA","ECOR66_PROKKA","ECOR6_PROKKA","ECOR8_PROKKA","KE21","KE24","KE25","KE26","KE4","KE40","KE41","KE46","KE47","KE48_PROKKA","KE50","KE51","KE54","KE55","KE58","LRPF007","LRPF62","MG1655","UTI89","W3110"]
    
    gene_id_to_info_map = {} 
    for index, row_data in df_roary.iterrows():
        roary_gene_name = row_data.get('Gene', f'UnnamedCluster_Row{index}')
        for strain_col_hdr in strain_names_list:
            if strain_col_hdr in row_data and pd.notna(row_data[strain_col_hdr]):
                for specific_gene_id_part in str(row_data[strain_col_hdr]).split(';'):
                    specific_gene_id = specific_gene_id_part.strip()
                    if specific_gene_id and specific_gene_id not in gene_id_to_info_map: 
                        gene_id_to_info_map[specific_gene_id] = {
                            'roary_gene_name': roary_gene_name,
                            'strain_name': strain_col_hdr
                        }
    
    all_cds_data = {}
    print("Pre-loading CDS files for sequence extraction...")
    strains_in_map = set(info['strain_name'] for info in gene_id_to_info_map.values())
    for s_name in strains_in_map: 
      cds_path = get_cds_file_path(s_name, args.cds_dir)
      cds_content_df = read_cds_file_to_df(cds_path)
      if cds_content_df is not None:
          all_cds_data[s_name] = cds_content_df
    print(f"Finished pre-loading {len(all_cds_data)} CDS files.")

    output_records = []
    for specific_gene_id, info_dict in gene_id_to_info_map.items():
        roary_name = info_dict['roary_gene_name']
        strain_identifier = info_dict['strain_name']
        sequence_data = "" 
        if strain_identifier in all_cds_data:
            cds_for_this_strain = all_cds_data[strain_identifier]
            seq_val = get_sequence_by_gene_id_from_df(specific_gene_id, cds_for_this_strain)
            if seq_val is not None:
                sequence_data = seq_val
        output_records.append([roary_name, specific_gene_id, strain_identifier, sequence_data])

    try:
        with open(args.output_csv, 'w', newline='', encoding='utf-8') as csvfile_out:
          writer_obj = csv.writer(csvfile_out)
          writer_obj.writerow(['gene', 'gene_id', 'strain_name', 'sequence'])
          writer_obj.writerows(output_records)
        print(f"Sequence data ({len(output_records)} records) written to {args.output_csv}")
    except Exception as e:
        print(f"Error writing output CSV {args.output_csv}: {e}")

if __name__ == '__main__':
    main()