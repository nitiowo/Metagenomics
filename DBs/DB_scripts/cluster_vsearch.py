#!/usr/bin/env python3
"""
cluster_by_tax_vsearch.py

Usage:
  python3 cluster_by_tax_vsearch.py --taxa_fasta taxa.fasta --id_fasta idseqs.fasta \
       --id_thresh 0.99 --out_prefix results --threads 2 [--keep_temp]

Produces:
 - results.cluster_report.tsv
 - results.taxonomy.fasta    (headers = taxonomy strings)
 - results.id.fasta          (headers = original ID-style headers)

Requires vsearch on PATH and Biopython installed.
"""
import argparse
import os
import sys
import tempfile
import subprocess
from collections import defaultdict, OrderedDict
from Bio import SeqIO

def check_vsearch():
    try:
        subprocess.run(["vsearch", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        sys.exit("ERROR: vsearch not found on PATH. Install vsearch and retry.")

def build_species_to_taxonomy_map(taxa_fasta):
    """Read taxonomy fasta where header is a semicolon-delimited taxonomy.
    Map species_name (last non-empty ; field) -> full taxonomy header (trimmed, with trailing ;).
    """
    sp2tax = {}
    for rec in SeqIO.parse(taxa_fasta, "fasta"):
        hdr = rec.description.strip()
        parts = [p.strip() for p in hdr.split(";") if p.strip() != ""]
        if not parts:
            continue
        species = parts[-1]
        taxstring = ";".join(parts) + ";"
        if species in sp2tax and sp2tax[species] != taxstring:
            # warn about conflicts (keep first)
            print(f"WARNING: species {species} has multiple taxonomy entries; keeping first.", file=sys.stderr)
            continue
        sp2tax[species] = taxstring
    return sp2tax

def extract_id_and_species_from_header(hdr):
    """From header like 'AADBE040-10 Glyptonotus antarcticus' return (id, species_str)"""
    parts = hdr.strip().split(None, 1)
    if len(parts) == 0:
        return ("", "")
    if len(parts) == 1:
        return (parts[0], "")
    return (parts[0], parts[1].strip())

def write_tmp_fasta_for_group(records_by_id, tmp_path):
    """Write a FASTA where each record header is the ID only (no spaces).
       records_by_id: OrderedDict id -> SeqRecord (original)
    """
    from Bio.SeqRecord import SeqRecord
    out_recs = []
    for rid, rec in records_by_id.items():
        # make a copy with id set to rid and empty description to ensure vsearch sees the id as header token
        r2 = SeqRecord(rec.seq, id=rid, description="")
        out_recs.append(r2)
    SeqIO.write(out_recs, tmp_path, "fasta")

def parse_uc_clusters(uc_path):
    """Parse vsearch .uc file to dict cluster_id -> list of member_ids.
       UC format columns: type, cluster, ... , query_label, target_label (tabs)
       We'll use col[0]=type, col[1]=cluster, col[8]=query_label, col[9]=target_label (if present).
    """
    clusters = defaultdict(list)
    # also track seeds (S lines)
    seeds = {}
    with open(uc_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            typ = cols[0]
            clid = cols[1]
            # query label in column 8 (0-based indexing)
            qlab = cols[8] if len(cols) > 8 else ""
            tlab = cols[9] if len(cols) > 9 else ""
            if typ == "S":
                # seed: qlab is the seed sequence name
                seeds[clid] = qlab
                clusters[clid].append(qlab)
            else:
                # H (hit) or other: qlab is member
                clusters[clid].append(qlab)
    return clusters, seeds

def choose_longest_member(cluster_members, id_to_seq):
    """Return the member id that has the longest sequence. id_to_seq maps id->sequence string."""
    best = None
    best_len = -1
    for m in cluster_members:
        seq = id_to_seq.get(m, "")
        l = len(seq)
        if l > best_len:
            best = m
            best_len = l
    return best, best_len

def main():
    p = argparse.ArgumentParser(description="Cluster sequences by taxonomy using vsearch and choose longest representative.")
    p.add_argument("--taxa_fasta", required=True, help="FASTA with taxonomy headers (semicolon-separated taxonomy in header).")
    p.add_argument("--id_fasta", required=True, help="FASTA with ID-style headers: >ID species_name")
    p.add_argument("--id_thresh", required=False, default=1.0, type=float, help="Identity threshold (e.g. 1.0 or 0.99).")
    p.add_argument("--out_prefix", required=True, help="Output prefix for report and fasta files.")
    p.add_argument("--threads", required=False, default=1, type=int, help="Threads passed to vsearch cluster (per-process).")
    p.add_argument("--keep_temp", action="store_true", help="Keep temporary files (for debugging).")
    args = p.parse_args()

    check_vsearch()

    taxa_fasta = args.taxa_fasta
    id_fasta = args.id_fasta
    id_thresh = float(args.id_thresh)
    out_prefix = args.out_prefix
    threads = int(args.threads)
    keep_temp = args.keep_temp

    print("Building species->taxonomy map from:", taxa_fasta, file=sys.stderr)
    sp2tax = build_species_to_taxonomy_map(taxa_fasta)
    print(f"Loaded {len(sp2tax)} species->taxonomy mappings.", file=sys.stderr)

    # Read id_fasta and group sequences by taxonomy
    print("Reading ID-style FASTA and grouping sequences by taxonomy...", file=sys.stderr)
    # We'll store for each taxonomy an OrderedDict of id->SeqRecord so order is stable
    tax_groups = defaultdict(OrderedDict)
    id_to_orig_header = {}
    id_to_seq = {}  # id -> sequence string
    unmapped_ids = []

    for rec in SeqIO.parse(id_fasta, "fasta"):
        orig_hdr = rec.description.strip()
        rid, species = extract_id_and_species_from_header(orig_hdr)
        if not rid:
            continue
        id_to_orig_header[rid] = orig_hdr
        seqstr = str(rec.seq).upper()
        id_to_seq[rid] = seqstr
        # map species -> taxonomy
        tax = sp2tax.get(species)
        if not tax:
            # try if the header itself is a taxonomy (semicolon present)
            if ";" in orig_hdr:
                # treat full header as taxonomy
                parts = [p.strip() for p in orig_hdr.split(";") if p.strip() != ""]
                if parts:
                    tax = ";".join(parts) + ";"
            # else unmapped
        if not tax:
            unmapped_ids.append(rid)
            continue
        # add record to that taxonomy group
        # store a SeqRecord with original sequence (we'll write tmp with ID-only header later)
        tax_groups[tax][rid] = rec

    print(f"Total tax groups found: {len(tax_groups)}", file=sys.stderr)
    if unmapped_ids:
        print(f"WARNING: {len(unmapped_ids)} sequences unmapped to taxonomy (they will be skipped).", file=sys.stderr)

    # Prepare outputs
    out_report = out_prefix + ".cluster_report.tsv"
    out_tax_fasta = out_prefix + ".taxonomy.fasta"
    out_id_fasta = out_prefix + ".id.fasta"

    rep_records_for_tax_fasta = []  # tuples (taxonomy_header, seq)
    rep_records_for_id_fasta = []   # tuples (orig_header, seq)

    # open report file
    rpt = open(out_report, "w")
    rpt.write("taxonomy\tcluster_id\trep_id\trep_len\tmember_count\tmember_ids\n")

    # process each taxonomy group
    group_count = len(tax_groups)
    gidx = 0
    for tax, recs_od in tax_groups.items():
        gidx += 1
        print(f"[{gidx}/{group_count}] Processing taxonomy: '{tax}'  (#seqs={len(recs_od)})", file=sys.stderr)
        if len(recs_od) == 0:
            continue

        # if only one sequence in group: choose it as representative with cluster_id 1
        if len(recs_od) == 1:
            rep_id = next(iter(recs_od.keys()))
            rep_seq = id_to_seq[rep_id]
            # write to outputs
            rep_records_for_tax_fasta.append((tax, rep_seq))
            rep_records_for_id_fasta.append((id_to_orig_header[rep_id], rep_seq))
            # report
            rpt.write(f"{tax}\t1\t{rep_id}\t{len(rep_seq)}\t1\t{rep_id}\n")
            continue

        # write temporary FASTA with ID-only headers (vsearch expects single-word sequence names)
        tmp_in = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
        tmp_in_name = tmp_in.name
        tmp_in.close()
        # write FASTA where header is just the ID (no spaces)
        write_tmp_fasta_for_group(recs_od, tmp_in_name)

        # sort by length
        tmp_sorted = tmp_in_name + ".sorted.fasta"
        cmd_sort = ["vsearch", "--sortbylength", tmp_in_name, "--output", tmp_sorted]
        subprocess.run(cmd_sort, check=True)

        # cluster
        tmp_centroids = tmp_in_name + ".centroids.fasta"
        tmp_uc = tmp_in_name + ".uc"
        cmd_cluster = [
            "vsearch", "--cluster_fast", tmp_sorted,
            "--id", str(id_thresh),
            "--centroids", tmp_centroids,
            "--uc", tmp_uc,
            "--threads", str(threads)
        ]
        subprocess.run(cmd_cluster, check=True)

        # parse uc to get clusters
        clusters, seeds = parse_uc_clusters(tmp_uc)
        # build mapping id->sequence (from recs_od)
        local_id_to_seq = {rid: id_to_seq[rid] for rid in recs_od.keys()}

        # For each cluster choose the absolute longest member
        # cluster IDs in uc are strings; sort cluster keys numerically if possible
        sorted_cluster_ids = sorted(clusters.keys(), key=lambda x: int(x) if x.isdigit() else x)

        for cid in sorted_cluster_ids:
            members = clusters[cid]
            # choose longest member
            chosen, chosen_len = choose_longest_member(members, local_id_to_seq)
            # append outputs
            rep_records_for_tax_fasta.append((tax, local_id_to_seq[chosen]))
            rep_records_for_id_fasta.append((id_to_orig_header[chosen], local_id_to_seq[chosen]))
            # write report line (member ids comma-separated)
            mlist = ",".join(members)
            rpt.write(f"{tax}\t{cid}\t{chosen}\t{chosen_len}\t{len(members)}\t{mlist}\n")

        # cleanup temporary files unless keep_temp
        if not keep_temp:
            for fn in (tmp_in_name, tmp_sorted, tmp_centroids, tmp_uc):
                try:
                    os.remove(fn)
                except OSError:
                    pass

    rpt.close()

    # write final FASTAs
    print("Writing output FASTAs...", file=sys.stderr)
    with open(out_tax_fasta, "w") as fh:
        for tax_hdr, seq in rep_records_for_tax_fasta:
            fh.write(f">{tax_hdr}\n")
            fh.write(seq + "\n")
    with open(out_id_fasta, "w") as fh:
        for orig_hdr, seq in rep_records_for_id_fasta:
            fh.write(f">{orig_hdr}\n")
            fh.write(seq + "\n")

    print("Done.")
    print("Report:", out_report)
    print("Taxonomy-style FASTA:", out_tax_fasta)
    print("ID-style FASTA:", out_id_fasta)

if __name__ == "__main__":
    main()
