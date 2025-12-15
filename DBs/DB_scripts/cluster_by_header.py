#!/usr/bin/env python3
"""
cluster_by_header.py

Cluster sequences within groups defined by their header format.

- Auto-detects header format (taxonomy-style if headers contain ';' or ID-style if 'ID species').
- Groups sequences by group_key:
    * taxonomy-style: group_key = canonical full taxonomy string (header normalized)
    * id-style:       group_key = species token (everything after first whitespace)
- For each group:
    * write temporary FASTA with ID-only headers
    * run vsearch --sortbylength (fallback to Python sort if fails)
    * run vsearch --cluster_fast --id <id_thresh> --centroids --uc
    * parse .uc, pick the absolute longest member per cluster
    * append representative sequences to output FASTA (header style matches input)
    * record a line in the cluster report TSV
- Cleans up temporary files by default (optionally keep with --keep_temp)

Outputs:
 - <out_prefix>.clustered.fasta  (representatives; header style matches input)
 - <out_prefix>.cluster_report.tsv (columns: group_key, cluster_id, rep_header, rep_len, member_count, member_ids_comma_sep)
 - optional: temporary files (if --keep_temp)

Requires:
  - Python3
  - biopython (pip install biopython)
  - vsearch on PATH (script uses fallback where possible)
"""

import argparse
import os
import sys
import tempfile
import subprocess
import shutil
import re
from collections import defaultdict, OrderedDict
from Bio import SeqIO

# ---------------------------
# Utility functions
# ---------------------------

def detect_header_format(sample_headers, verbose=False):
    """
    Decide whether headers look like taxonomy-style or id-style.
    Returns 'taxonomy' or 'id'.
    Heuristic:
      - if many headers contain ';' -> taxonomy
      - else if many headers have a first token that looks like an ID and content after a space -> id
      - fallback to taxonomy if semicolons present else id
    """
    n = len(sample_headers)
    if n == 0:
        return "id"
    semis = sum(1 for h in sample_headers if ';' in h)
    spaces = sum(1 for h in sample_headers if re.search(r'\s+', h))
    if verbose:
        print(f"[detect] sampled {n} headers: semicolons={semis}, spaces={spaces}", file=sys.stderr)
    # If at least ~60% have semicolons prefer taxonomy
    if semis >= max(1, int(0.6 * n)):
        return "taxonomy"
    # else if many have spaces and few semicolons, choose id
    if spaces >= max(1, int(0.6 * n)) and semis < (0.3 * n):
        return "id"
    # fallback: if any semicolon present -> taxonomy else id
    return "taxonomy" if semis > 0 else "id"

def canonicalize_taxonomy_header(hdr):
    """
    Convert a taxonomy header into a canonical form:
      - split on ';' or ','; strip fields; rejoin with ';' and ensure trailing ';'
    """
    parts = [p.strip() for p in re.split(r'[;,]', hdr) if p.strip() != ""]
    if not parts:
        return hdr.strip()
    return ";".join(parts) + ";"

def parse_header_for_grouping(header, mode):
    """
    Parse header string (no leading '>').
    Returns tuple (record_id_token, group_key, canonical_output_header).

    record_id_token: single-word id used in temporary FASTA (unique if duplicates exist)
    group_key: grouping key used when aggregating records (taxonomy string or species token)
    canonical_output_header: header string to use when writing representative output
    """
    hdr = header.strip()
    if mode == "taxonomy":
        # treat the whole header as taxonomy (canonicalized)
        tax = canonicalize_taxonomy_header(hdr)
        # attempt to make a unique record id (take last token as id-like if present, else use
        # a short hash). But better to use the first token (up to first space) if that exists.
        first_token = hdr.split(None, 1)[0]
        rec_id = first_token
        # ensure rec_id is safe (no spaces)
        rec_id = re.sub(r'\s+', '_', rec_id)
        # canonical header (taxonomy-style) is tax
        return (rec_id, tax, tax)
    else:
        # ID-style header like "AADBE040-10 Glyptonotus antarcticus"
        parts = hdr.split(None, 1)
        rec_id = parts[0] if parts else ""
        species = parts[1].strip() if len(parts) > 1 else ""
        # canonical species key: normalize whitespace
        species_key = re.sub(r'\s+', ' ', species).strip()
        # canonical output header: keep original header (so we can write same format)
        return (rec_id, species_key, hdr)

def write_temp_fasta_idonly(records_od, tmp_path):
    """
    Write a temporary FASTA with headers equal to the record ids (single-token)
    records_od: OrderedDict mapping record_id -> SeqRecord (original SeqRecord)
    """
    # Create SeqRecord objects with id = record_id and description empty (vsearch expects single-word id)
    out_recs = []
    from Bio.SeqRecord import SeqRecord
    for rid, rec in records_od.items():
        # ensure sequence is a string; create SeqRecord with id only
        r2 = SeqRecord(rec.seq, id=rid, description="")
        out_recs.append(r2)
    SeqIO.write(out_recs, tmp_path, "fasta")

def python_sortbylength(in_fasta, out_fasta):
    """Fallback sorting by sequence length in Python (descending)."""
    recs = list(SeqIO.parse(in_fasta, "fasta"))
    recs.sort(key=lambda r: len(r.seq), reverse=True)
    SeqIO.write(recs, out_fasta, "fasta")

def run_vsearch_sort(in_fasta, out_sorted, tmpdir=None):
    """Try vsearch --sortbylength, return True on success, False on failure."""
    cmd = ["vsearch", "--sortbylength", in_fasta, "--output", out_sorted]
    if tmpdir:
        # vsearch supports --tmpdir on some builds; but not universal. Instead set environment TMPDIR.
        env = os.environ.copy()
        env["TMPDIR"] = tmpdir
    else:
        env = None
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        return True
    except Exception as e:
        # vsearch sort failed; caller should fallback to python sort
        return False

def run_vsearch_cluster(sorted_fasta, centroids_path, uc_path, id_thresh, threads=1):
    """
    Run vsearch clustering. Raises CalledProcessError if vsearch fails.
    """
    cmd = [
        "vsearch", "--cluster_fast", sorted_fasta,
        "--id", str(id_thresh),
        "--centroids", centroids_path,
        "--uc", uc_path,
        "--threads", str(threads)
    ]
    subprocess.run(cmd, check=True)

def parse_uc_file(uc_path):
    """
    Parse vsearch .uc file.
    Returns: clusters dict cluster_id -> list(member_labels)
             seeds dict cluster_id -> seed_label (if present)
    UC columns: type, cluster, size, perc, strand, qstart, qend, qlen, qlabel, tlabel (tab-delimited)
    We use qlabel (col index 8)
    """
    clusters = defaultdict(list)
    seeds = {}
    with open(uc_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            typ = cols[0]
            clid = cols[1]
            # qlabel is usually column 8 (0-based idx)
            qlabel = cols[8] if len(cols) > 8 else ""
            if typ == "S":
                # S line is a seed (cluster representative)
                seeds[clid] = qlabel
                clusters[clid].append(qlabel)
            else:
                # H or other -> a member
                clusters[clid].append(qlabel)
    return clusters, seeds

def choose_longest_member(members, id_to_seq):
    """
    Given list of member ids and id->sequence mapping, choose the member with the
    longest sequence. Returns (chosen_id, length).
    """
    best = None
    bestlen = -1
    for m in members:
        seq = id_to_seq.get(m, "")
        L = len(seq)
        if L > bestlen:
            best = m
            bestlen = L
    return best, bestlen

# ---------------------------
# Main clustering flow
# ---------------------------

def cluster_by_header(in_fasta, id_thresh, out_prefix, threads=1, tmpdir=None, keep_temp=False, wrap=None, sample_headers=100):
    """
    Main coordinator function.
    - Detect header style
    - Group records by group_key
    - For each group: cluster with vsearch, pick longest reps, write outputs & report
    """
    # 1) sample headers to detect format
    sample = []
    with open(in_fasta) as fh:
        # read up to sample_headers records cheaply by scanning lines
        for line in fh:
            if line.startswith(">"):
                sample.append(line[1:].strip())
                if len(sample) >= sample_headers:
                    break
    mode = detect_header_format(sample, verbose=True)
    print(f"[info] detected header mode = {mode}", file=sys.stderr)

    # 2) group records by group_key
    print("[info] grouping records by group_key ...", file=sys.stderr)
    groups = defaultdict(OrderedDict)   # group_key -> OrderedDict(record_id -> SeqRecord)
    id_to_seq = {}                      # record_id -> sequence string
    id_to_orig_header = {}              # record_id -> original header string (no '>')
    seen_ids = defaultdict(int)         # help ensure unique rec ids if duplicates exist
    unmapped = []

    # iterate through FASTA streaming
    for rec in SeqIO.parse(in_fasta, "fasta"):
        orig_hdr = rec.description.strip()
        rec_id_token, group_key, canonical_hdr = parse_header_for_grouping(orig_hdr, mode)
        # ensure a unique rec_id (if duplicate rec_id_token appear)
        count = seen_ids[rec_id_token]
        seen_ids[rec_id_token] += 1
        if count > 0:
            rec_id = f"{rec_id_token}_{count}"
        else:
            rec_id = rec_id_token if rec_id_token else f"rec_{len(id_to_seq)+1}"
        # store
        groups[group_key][rec_id] = rec
        id_to_seq[rec_id] = str(rec.seq).upper()
        id_to_orig_header[rec_id] = orig_hdr

    total_groups = len(groups)
    print(f"[info] grouped into {total_groups} groups", file=sys.stderr)

    # Prepare outputs: open files for streaming
    out_fasta = out_prefix + ".clustered.fasta"
    report_tsv = out_prefix + ".cluster_report.tsv"
    unmapped_fasta = out_prefix + ".unmapped.fasta"

    out_fh = open(out_fasta, "w")
    rpt_fh = open(report_tsv, "w")
    rpt_fh.write("group_key\tcluster_id\trep_header\trep_len\tmember_count\tmember_ids\n")

    # process each group sequentially
    gcount = 0
    total_reps = 0
    for gkey, recs_od in groups.items():
        gcount += 1
        nseqs = len(recs_od)
        print(f"[{gcount}/{total_groups}] Processing group: '{gkey}'  (#seqs={nseqs})", file=sys.stderr)

        # if group has 0 sequences skip
        if nseqs == 0:
            continue

        # If only one sequence, choose it as representative without running vsearch
        if nseqs == 1:
            rep_id = next(iter(recs_od.keys()))
            rep_seq = id_to_seq[rep_id]
            # determine header to write based on mode:
            if mode == "taxonomy":
                # header is the group taxonomy string (gkey)
                out_hdr = gkey
            else:
                # header is the original header for the rep
                out_hdr = id_to_orig_header[rep_id]
            # write representative to output fasta (sequence on single line)
            out_fh.write(f">{out_hdr}\n{rep_seq}\n")
            # write to report
            rpt_fh.write(f"{gkey}\t1\t{out_hdr}\t{len(rep_seq)}\t1\t{rep_id}\n")
            total_reps += 1
            continue

        # For groups with multiple sequences: prepare temporary FASTA
        tmp_in = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", dir=tmpdir)
        tmp_in_name = tmp_in.name
        tmp_in.close()
        try:
            # write temp fasta with ID-only headers
            write_temp_fasta_idonly(recs_od, tmp_in_name)

            # sorted output temp path
            tmp_sorted = tmp_in_name + ".sorted.fasta"
            # try vsearch sort, fallback to python sort
            ok = False
            try:
                ok = run_vsearch_sort(tmp_in_name, tmp_sorted, tmpdir=tmpdir)
            except Exception:
                ok = False
            if not ok:
                # fallback
                print(f"[warn] vsearch --sortbylength failed for group; using Python fallback sort", file=sys.stderr)
                python_sortbylength(tmp_in_name, tmp_sorted)

            # run clustering with vsearch
            tmp_centroids = tmp_in_name + ".centroids.fasta"
            tmp_uc = tmp_in_name + ".uc"
            try:
                run_vsearch_cluster(tmp_sorted, tmp_centroids, tmp_uc, id_thresh, threads=threads)
            except subprocess.CalledProcessError as e:
                # vsearch cluster failed; fallback to exact dedup by sequence string
                print(f"[error] vsearch --cluster_fast failed for group {gkey}: {e}. Falling back to exact dedup.", file=sys.stderr)
                # implement exact dedup fallback: group identical sequences and pick one (longest)
                seq_to_ids = defaultdict(list)
                for rid in recs_od.keys():
                    seq_to_ids[id_to_seq[rid]].append(rid)
                # for each unique sequence choose one representative (longest is same here), append to outputs and report
                cid = 0
                for seq, members in seq_to_ids.items():
                    cid += 1
                    chosen = members[0]
                    # header as above
                    out_hdr = gkey if mode == "taxonomy" else id_to_orig_header[chosen]
                    out_fh.write(f">{out_hdr}\n{seq}\n")
                    rpt_fh.write(f"{gkey}\t{cid}\t{out_hdr}\t{len(seq)}\t{len(members)}\t{','.join(members)}\n")
                    total_reps += 1
                # cleanup tmp files later and continue to next group
                if not keep_temp:
                    try:
                        os.remove(tmp_in_name)
                        os.remove(tmp_sorted)
                    except OSError:
                        pass
                continue

            # parse .uc to get clusters
            clusters, seeds = parse_uc_file(tmp_uc)
            # For each cluster choose longest member
            # sort cluster ids for deterministic order
            try:
                cluster_keys = sorted(clusters.keys(), key=lambda x: int(x))
            except Exception:
                cluster_keys = sorted(clusters.keys())

            for cid in cluster_keys:
                members = clusters[cid]
                chosen, chosenlen = choose_longest_member(members, id_to_seq)
                if not chosen:
                    # defensive: skip empty cluster
                    continue
                # determine header to write
                if mode == "taxonomy":
                    out_hdr = gkey
                else:
                    out_hdr = id_to_orig_header.get(chosen, chosen)
                out_seq = id_to_seq.get(chosen, "")
                # write representative (single-line seq)
                out_fh.write(f">{out_hdr}\n{out_seq}\n")
                # write report row
                rpt_fh.write(f"{gkey}\t{cid}\t{out_hdr}\t{chosenlen}\t{len(members)}\t{','.join(members)}\n")
                total_reps += 1

        finally:
            # cleanup temp files unless keep_temp
            if keep_temp:
                print(f"[debug] kept temp files for group {gkey}: {tmp_in_name}*", file=sys.stderr)
            else:
                for fn in (tmp_in_name, tmp_sorted, tmp_centroids, tmp_uc):
                    try:
                        if os.path.exists(fn):
                            os.remove(fn)
                    except OSError:
                        pass

    out_fh.close()
    rpt_fh.close()

    print(f"[done] processed {gcount} groups; wrote {total_reps} representative sequences", file=sys.stderr)
    print(f"[done] fasta -> {out_fasta}", file=sys.stderr)
    print(f"[done] report -> {report_tsv}", file=sys.stderr)
    if unmapped:
        with open(unmapped_fasta, "w") as ufh:
            for rid in unmapped:
                seq = id_to_seq.get(rid, "")
                hdr = id_to_orig_header.get(rid, rid)
                ufh.write(f">{hdr}\n{seq}\n")
        print(f"[info] some records were unmapped; see {unmapped_fasta}", file=sys.stderr)

# ---------------------------
# CLI
# ---------------------------

def cli():
    p = argparse.ArgumentParser(description="Cluster sequences within header-defined groups (taxonomy-style or ID-style headers).")
    p.add_argument("--in_fasta", required=True, help="Input FASTA (taxonomy-style or ID-style headers).")
    p.add_argument("--id_thresh", type=float, default=1.0, help="Identity threshold for vsearch clustering (e.g. 1.0, 0.99).")
    p.add_argument("--out_prefix", required=True, help="Output prefix for clustered fasta and report.")
    p.add_argument("--threads", type=int, default=1, help="Threads to pass to vsearch.")
    p.add_argument("--tmpdir", default=None, help="Directory to store temporary files (default: system temp).")
    p.add_argument("--keep_temp", action="store_true", help="Keep temporary files for debugging (default: delete).")
    p.add_argument("--wrap", type=int, default=0, help="Wrap sequence lines to this width (0 = no wrap, default).")
    args = p.parse_args()

    # basic checks
    if not os.path.exists(args.in_fasta):
        sys.exit(f"ERROR: input fasta not found: {args.in_fasta}")
    # Check vsearch exists (we will still fallback sort if sort fails, but clustering requires vsearch)
    try:
        subprocess.run(["vsearch", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("[warn] vsearch not found or not runnable; clustering step will likely fail unless vsearch is available.", file=sys.stderr)

    cluster_by_header(
        in_fasta=args.in_fasta,
        id_thresh=args.id_thresh,
        out_prefix=args.out_prefix,
        threads=args.threads,
        tmpdir=args.tmpdir,
        keep_temp=args.keep_temp,
        wrap=args.wrap
    )

if __name__ == "__main__":
    cli()
