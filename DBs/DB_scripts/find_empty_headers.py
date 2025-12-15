#!/usr/bin/env python
from Bio import SeqIO
import sys

def find_empty_fasta_records(filename):
    """
    Parses a FASTA file and prints the headers of empty records.
    """
    empty_count = 0
    print(f"Checking file: {filename}")
    
    try:
        # Use SeqIO.parse to iterate through records
        for record in SeqIO.parse(filename, "fasta"):
            if len(record.seq) == 0:
                empty_count += 1
                # Output the full header/description of the empty record
                print(f"Empty record found: {record.description}")
                # You can also use record.id for just the part before the first space
                # print(f"Empty record ID: {record.id}")

        print(f"\nTotal number of empty records found: {empty_count}")

    except IOError as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python find_empty_headers.py <your_file.fasta>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    find_empty_fasta_records(fasta_file)
