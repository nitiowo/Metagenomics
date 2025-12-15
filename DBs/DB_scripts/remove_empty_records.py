#!/usr/bin/env python
from Bio import SeqIO
import sys
import os

def remove_empty_records(input_filename, output_filename):
    # ... (function body remains the same as before) ...
    """
    Reads a FASTA file, removes records with empty sequences,
    and writes the remaining records to a new file.
    """
    count_before = 0
    count_after = 0
    non_empty_records = []

    try:
        # Pass the filename string to SeqIO.parse
        for record in SeqIO.parse(input_filename, "fasta"):
            count_before += 1
            if len(record.seq) > 0:
                non_empty_records.append(record)
                count_after += 1
        
        with open(output_filename, "w") as output_handle:
            SeqIO.write(non_empty_records, output_handle, "fasta")
            
        print(f"Processed file: {input_filename}")
        print(f"Total records found: {count_before}")
        print(f"Records removed (empty): {count_before - count_after}")
        print(f"Records retained (non-empty): {count_after}")
        print(f"Filtered file saved to: {output_filename}")

    except IOError as e:
        print(f"Error handling file(s): {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_empty.py <input_file.fasta> <output_file.fasta>")
        sys.exit(1)
    
    # Corrected lines: Access the specific file path strings from the list
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    remove_empty_records(input_file_path, output_file_path)
