#!/usr/bin/env python3
# filter_relabel_fasta.py
# Usage: python filter_relabel_fasta.py refs.csv input.fasta output.fasta

import sys
import csv
from Bio import SeqIO

csv_path = sys.argv[1]
fasta_in  = sys.argv[2]
fasta_out = sys.argv[3]

# Build mapping ID -> LevelHeader from CSV (streamed, minimal memory)
# Adjust column names below if your CSV headers differ
id_col = "ID"
cols = ["Kingdom","Phylum","Class","Order","Family","Genus"]  # levels 1..6

map_header = {}
with open(csv_path, newline='') as fh:
    rdr = csv.DictReader(fh)
    for row in rdr:
        rid = row.get(id_col)
        if not rid:
            continue
        # build Level1;Level2;...;Level6; string (coalesce missing to empty)
        header = ";".join((row.get(c,"").strip() for c in cols)) + ";"
        map_header[rid.strip()] = header

# Stream FASTA and write matches
with open(fasta_out, "w") as outfh:
    for rec in SeqIO.parse(fasta_in, "fasta"):
        # rec.description contains the full header line (without leading '>')
        # the first '|' field is the ID
        fasta_id = rec.description.split("|", 1)[0]
        if fasta_id in map_header:
            rec.id = map_header[fasta_id]    # SeqIO will write >rec.id
            rec.description = ""             # avoid adding original desc
            SeqIO.write(rec, outfh, "fasta")