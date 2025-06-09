# blast_central_genes.py
import ast
import os
import csv
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse

def extract_record_from_cds_file(cds_filepath_val, target_locus_id_val, strain_name_val):
    """
    Extracts a FASTA record for a specific locus ID from a given CDS file.
    CDS file is expected to be tab-delimited:
    col 0: accession, col 1: start, col 2: end, col 3: locus_id/gene_id,
    col 4: N/A, col 5: strand, col 6: gene_name/description, col 7: sequence
    """
    try:
        with open(cds_filepath_val, 'r', encoding='utf-8') as f_cds: # 
            for line_cds in f_cds: # 
                line_cds = line_cds.strip() # 
                if not line_cds: continue # 
                fields_cds = line_cds.split("\t") # 
                if len(fields_cds) < 8: continue # Ensure enough columns 
                # Compare target_locus_id_val (from reversed_df) with 4th column (index 3) of CDS file
                if fields_cds[3].strip() == str(target_locus_id_val): # 
                    header_val = f"{strain_name_val}|{fields_cds[6].strip()}|{target_locus_id_val}" # Use 7th col (idx 6) for description 
                    sequence_val = fields_cds[7].strip() # 8th col (idx 7) is sequence # 
                    return SeqRecord(Seq(sequence_val), id=header_val, description="") # 
        return None # Locus ID not found 
    except FileNotFoundError: # 
        # print(f"Warning: CDS file not found {cds_filepath_val}")
        return None # 
    except Exception as e_cds: # 
        print(f"Error processing CDS file {cds_filepath_val} for locus {target_locus_id_val}: {e_cds}") # 
        return None # 

def main(): # 
    parser = argparse.ArgumentParser(description="BLASTs central genes against strains where they are reported missing.") # 
    parser.add_argument('--input_reversed_csv', required=True, help="Path to the input CSV (e.g., reversed_df_completeid.csv) which contains 'central_gene', 'strain', and 'position_0' columns with tuple data.") # 
    parser.add_argument('--cds_files_dir', required=True, help="Directory containing CDS sequence files (e.g., <strain_name>.CDS.txt).") # 
    parser.add_argument('--genome_fasta_dir', required=True, help="Directory containing whole genome FASTA files (e.g., <strain_name>.fasta).") # 
    parser.add_argument('--blast_output_dir', required=True, help="Directory to save all BLAST-related output files (queries, DBs, results).") # 
    parser.add_argument('--no_hits_output_file', required=True, help="Path to the output CSV file for listing gene-strain pairs with no BLAST hits.") # 
    parser.add_argument('--blast_bin_path', help="Optional: Full path to the BLAST+ 'bin' directory if not in system PATH.") # 
    args = parser.parse_args() # 

    # Update PATH for BLAST+ if a specific path is provided
    if args.blast_bin_path and args.blast_bin_path.strip() != "": # 
        blast_exec_path = os.path.expanduser(args.blast_bin_path) # 
        current_path = os.environ.get("PATH", "") # 
        os.environ["PATH"] = blast_exec_path + os.pathsep + current_path # 
        print(f"Updated system PATH to include BLAST+ binaries from: {blast_exec_path}") # 

    # Ensure the main BLAST output directory exists
    final_results_dir = args.blast_output_dir # 
    if not os.path.isdir(final_results_dir): # 
        try: # 
            os.makedirs(final_results_dir) # 
            print(f"Created BLAST output directory: {final_results_dir}") # 
        except OSError as e: # 
            print(f"Error: Could not create BLAST output directory {final_results_dir}: {e}") # 
            return # 

    # Process the input CSV to group genes and identify present/missing strains
    gene_info_groups = {} # 
    try: # 
        with open(args.input_reversed_csv, 'r', encoding='utf-8') as f_csv_in: # 
            reader = csv.DictReader(f_csv_in) # 
            for r_csv in reader: # 
                gene_name_key = r_csv.get("central_gene", "").strip().lower() # 
                strain_name_val = r_csv.get("strain", "").strip() # 
                if not gene_name_key or not strain_name_val: continue # Skip rows with missing essential info # 

                if gene_name_key not in gene_info_groups: # 
                    gene_info_groups[gene_name_key] = {"present": [], "missing": []} # 
                
                position_0_val = r_csv.get("position_0", "").strip()
                # Modified condition to better detect missing genes
                if (not position_0_val or 
                    position_0_val.lower() == "nan" or 
                    position_0_val == "" or 
                    position_0_val == "[]" or 
                    position_0_val == "()" or
                    position_0_val == "{}"):
                    gene_info_groups[gene_name_key]["missing"].append(strain_name_val)
                    print(f"Found missing gene {gene_name_key} in strain {strain_name_val}")
                else:
                    gene_info_groups[gene_name_key]["present"].append((r_csv, strain_name_val)) # 
        # print(f"Processed input CSV: {args.input_reversed_csv}")
    except FileNotFoundError: # 
        print(f"Error: Input CSV file '{args.input_reversed_csv}' not found.") # 
        return # 
    except Exception as e_csv_read: # 
        print(f"Error reading CSV file '{args.input_reversed_csv}': {e_csv_read}") # 
        return # 

    print("\nDEBUG: Input Data Analysis")
    print(f"Total rows read: {sum(1 for _ in csv.DictReader(open(args.input_reversed_csv)))}")
    print("\nSample position_0 values:")
    with open(args.input_reversed_csv, 'r') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i < 5:  # Print first 5 rows
                print(f"Gene: {row['central_gene']}, Strain: {row['strain']}, Position_0: {row['position_0']}")

    query_fasta_files_map = {} # Maps central_gene_name to its query FASTA file path # 
    missing_strains_map = {}   # Maps central_gene_name to list of strains where it's "missing" # 

    print("--- Generating Query FASTA files ---")
    for gene_key_name, g_data in gene_info_groups.items(): # 
        if not g_data["present"]: # 
            # print(f"DEBUG: No 'present' entries for gene '{gene_key_name}' to create query FASTA.")
            continue # 
        # print(f"Processing gene: {gene_key_name} with {len(g_data['present'])} present strains") #
        query_file_path = os.path.join(final_results_dir, f"query_{gene_key_name}.fasta") # 
        fasta_records_for_query = [] # 

        for (csv_row_data, strain_val) in g_data["present"]: # 
            extracted_locus_id = None # 
            pos0_cell_content = csv_row_data.get("position_0", "") # This column should contain the tuple string # 

            if isinstance(pos0_cell_content, str) and pos0_cell_content.strip() and pos0_cell_content.lower() != "nan": # 
                try: # 
                    # The tuple from create_sorted_strains_table is:
                    # (GFF_gene_name, GFF_ID, strand, core_status, start, end)
                    parsed_tuple = ast.literal_eval(pos0_cell_content) # 
                    if isinstance(parsed_tuple, tuple) and len(parsed_tuple) >= 2: # Need at least GFF_gene and GFF_ID # 
                         extracted_locus_id = str(parsed_tuple[1]).strip() # GFF ID is at index 1 of the tuple # 
                except (ValueError, SyntaxError): # 
                    # print(f"Warning: Could not parse position_0 content as tuple for {gene_key_name} in {strain_val}: '{pos0_cell_content}'")
                    pass # 
            
            if extracted_locus_id is None: # 
                # print(f"Warning: No valid locus ID extracted from position_0 for {gene_key_name} in strain {strain_val}. Skipping this record for query FASTA.")
                continue # 
            
            cds_f_path = os.path.join(args.cds_files_dir, f"{strain_val}.CDS.txt") # 
            if not os.path.isfile(cds_f_path): # 
                # print(f"Warning: CDS file not found for strain {strain_val}: {cds_f_path}")
                continue # 
            
            seq_record = extract_record_from_cds_file(cds_f_path, extracted_locus_id, strain_val) # 
            if seq_record: # 
                fasta_records_for_query.append(seq_record) # 
            # else:
                # print(f"Warning: Sequence record not extracted for locus {extracted_locus_id} from {cds_f_path}") # 
        
        if fasta_records_for_query: # 
            try: # 
                with open(query_file_path, 'w', encoding='utf-8') as f_query_out: # 
                    SeqIO.write(fasta_records_for_query, f_query_out, "fasta") # 
                query_fasta_files_map[gene_key_name] = query_file_path # 
                print(f"DEBUG: Writing {len(fasta_records_for_query)} sequences to query FASTA: {query_file_path}") # 
            except IOError as e_io_q: # 
                 print(f"Error writing query FASTA for '{gene_key_name}' to '{query_file_path}': {e_io_q}") # 
                 continue # Skip this gene if query FASTA cannot be written # 
        # else:
            # print(f"DEBUG: No sequences collected for query FASTA for gene '{gene_key_name}'.") # 

        missing_strains_map[gene_key_name] = g_data["missing"] # 
        print(f"DEBUG: Gene '{gene_key_name}' has {len(g_data['missing'])} missing strains.") #

    no_blast_hits_found_list = [] # 
    if not missing_strains_map or all(not v for v in missing_strains_map.values()): # Check if map is empty or all lists of missing strains are empty
        print("--- No missing strains identified for any gene with a query FASTA. BLAST searches will be effectively skipped. ---")
    else:
        print(f"--- Starting BLAST Operations ({len(missing_strains_map)} gene groups to process) ---") # 

    for gene_blast_key, strains_to_blast_list in missing_strains_map.items(): # 
        current_query_fasta = query_fasta_files_map.get(gene_blast_key) # 
        if not current_query_fasta: # 
            # print(f"DEBUG: No query FASTA for gene '{gene_blast_key}', skipping BLAST.")
            continue # 
        
        if not strains_to_blast_list: # If the list of missing strains for this gene is actually empty
            # print(f"DEBUG: No missing strains listed for gene '{gene_blast_key}' in the input CSV logic.")
            continue

        print(f"Processing gene: {gene_blast_key} (Query: {current_query_fasta})") # 
        for missing_strain_name in strains_to_blast_list: # 
            print(f"  Targeting missing strain: {missing_strain_name}") # 
            target_genome_fasta = os.path.join(args.genome_fasta_dir, f"{missing_strain_name}.fasta") # 
            if not os.path.isfile(target_genome_fasta): # 
                print(f"  Warning: Genome FASTA for {missing_strain_name} not found at {target_genome_fasta}") # 
                continue # 

            blast_db_prefix = os.path.join(final_results_dir, f"db_{missing_strain_name}") # 
            makeblastdb_log_file = os.path.join(final_results_dir, f"makeblastdb_{missing_strain_name}.log") # 
            
            db_exists_check_file = blast_db_prefix + ".nhr" # 
            if not os.path.exists(db_exists_check_file): # 
                makeblastdb_cmd = ["makeblastdb", "-in", target_genome_fasta, "-dbtype", "nucl", "-out", blast_db_prefix, "-logfile", makeblastdb_log_file] # 
                print(f"    Creating BLAST DB for {missing_strain_name}...") # 
                try: # 
                    subprocess.run(makeblastdb_cmd, check=True, capture_output=True, text=True, encoding='utf-8') # 
                except subprocess.CalledProcessError as e_db: # 
                    print(f"    ERROR: makeblastdb failed for {missing_strain_name}. Command: {' '.join(e_db.cmd)}. Stderr: {e_db.stderr.strip()}") # 
                    continue # Skip to next strain if DB creation fails # 
            # else:
                # print(f"    Info: BLAST DB for {missing_strain_name} already exists.") # 

            # Verbose BLAST output
            blast_verbose_out_path = os.path.join(final_results_dir, f"blast_{gene_blast_key}_vs_{missing_strain_name}.verbose.out") # 
            blastn_verbose_cmd = ["blastn", "-query", current_query_fasta, "-db", blast_db_prefix, "-out", blast_verbose_out_path] # 
            print(f"    Running BLASTN (verbose) for {gene_blast_key} vs {missing_strain_name}...") # 
            try: # 
                subprocess.run(blastn_verbose_cmd, check=True, capture_output=True, text=True, encoding='utf-8') # 
            except subprocess.CalledProcessError as e_blastn_v: # 
                print(f"    ERROR: Verbose BLASTN failed for {gene_blast_key} vs {missing_strain_name}. Stderr: {e_blastn_v.stderr.strip()}") # 
                # Continue, but tabular might also fail or be meaningless
            
            # Check verbose output for "No hits found"
            hit_actually_found_in_verbose = False # Assume no hit initially
            try: # 
                with open(blast_verbose_out_path, 'r', encoding='utf-8') as bf_read: # 
                    if "***** No hits found *****" not in bf_read.read(): # If "No hits found" is NOT present, then hits WERE found
                        hit_actually_found_in_verbose = True
                    else: # 
                        no_blast_hits_found_list.append((gene_blast_key, missing_strain_name)) # 
                        print(f"    Info: No BLAST hits found for gene '{gene_blast_key}' in strain '{missing_strain_name}' (verbose).") # 
            except FileNotFoundError: # 
                 print(f"    Warning: BLAST verbose output file not found: {blast_verbose_out_path}. Assuming no hits or prior error.") # 
                 if not (gene_blast_key, missing_strain_name) in no_blast_hits_found_list: # Avoid adding duplicate no-hit entries
                    no_blast_hits_found_list.append((gene_blast_key, missing_strain_name)) # 
            except Exception as e_read_blast_v: # 
                print(f"    Error reading {blast_verbose_out_path}: {e_read_blast_v}") # 
                # Treat as no verifiable hit if read error occurs
                if not (gene_blast_key, missing_strain_name) in no_blast_hits_found_list:
                    no_blast_hits_found_list.append((gene_blast_key, missing_strain_name))


            # Tabular BLAST output
            if hit_actually_found_in_verbose : # 
                blast_tabular_out_path = os.path.join(final_results_dir, f"blast_{gene_blast_key}_vs_{missing_strain_name}.tabular.tsv") # 
                outfmt_string = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle" # 
                blastn_tabular_cmd = [ # 
                    "blastn", "-query", current_query_fasta, "-db", blast_db_prefix, # 
                    "-max_target_seqs", "1", "-max_hsps", "1",  # 
                    "-outfmt", outfmt_string, "-out", blast_tabular_out_path # 
                ] # 
                print(f"    Running BLASTN (tabular) for {gene_blast_key} vs {missing_strain_name}...") # 
                try: # 
                    subprocess.run(blastn_tabular_cmd, check=True, capture_output=True, text=True, encoding='utf-8') # 
                    
                    header_line = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tslen\tstitle\n" # 
                    try: # 
                        with open(blast_tabular_out_path, 'r+', encoding='utf-8') as f_in_tab: # r+ to read and then overwrite # 
                            content_tab = f_in_tab.read() # 
                            f_in_tab.seek(0, 0) # Go to the beginning of the file # 
                            f_in_tab.write(header_line) # 
                            if content_tab.strip(): # Write content only if it's not empty # 
                                 f_in_tab.write(content_tab) # 
                    except FileNotFoundError: pass  # 
                    except IOError: print(f"    IOError with {blast_tabular_out}") # 
                except subprocess.CalledProcessError as e_blastn_t: # 
                    print(f"    ERROR: Tabular BLASTN failed for {gene_blast_key} vs {missing_strain_name}. Stderr: {e_blastn_t.stderr.strip()}") # 
            
            # Cleanup BLAST database files for this strain 
            # Initialize db_files_to_try_remove before using it
            db_files_to_try_remove = []  # Add this line before appending to the list

            # Cleanup BLAST database files for this strain
            db_extensions_to_remove = ['.nhr', '.nin', '.nsq', '.ndb', '.nos', '.njs', '.nal', '.nni', '.nnd', '.log']  # Extensions to remove
            db_files_to_try_remove.append(blast_db_prefix)  # Append the prefix to the list

            for db_file in db_files_to_try_remove:  # Iterate over files to remove
                for ext in db_extensions_to_remove:  # Iterate over extensions
                    file_to_remove = db_file + ext  # Construct the full file path
                    if os.path.exists(file_to_remove):  # Check if the file exists
                        try:
                            os.remove(file_to_remove)  # Attempt to remove the file
                            print(f"Removed BLAST DB file: {file_to_remove}")
                        except OSError as e:
                            print(f"Error removing BLAST DB file {file_to_remove}: {e}")

    print("Finished all BLAST operations.") # 

    try: # 
        with open(args.no_hits_output_file, 'w', newline='', encoding='utf-8') as f_no_hits: # 
            writer = csv.writer(f_no_hits) # 
            writer.writerow(['Gene', 'Strain']) # Write header # 
            if no_blast_hits_found_list: # 
                writer.writerows(no_blast_hits_found_list) # 
        print(f"\nList of gene-strain pairs with no BLAST hits saved to: {args.no_hits_output_file}") # 
    except IOError as e_io_nh: # 
        print(f"\nError writing no-hits list to file {args.no_hits_output_file}: {e_io_nh}") # 

    print(f"BLAST outputs are in: {final_results_dir}") # 
    print("BLAST script completed.") # 

if __name__ == '__main__': # 
    main() #