#!/usr/bin/env python3
"""Renames the cluster representatives' FASTA to <lineage ID>.<n>, the lineage ID being the cluster's LCA taxid.
Inputs: the LCA DB (mmseqs lca output; entries "taxid\trank\tname\n\0", data split in parts .0 .1 ..., offsets global), the
representatives' sub-DB (its index gives the keys in the order mmseqs convert2fasta wrote them) and that FASTA. The two indexes
are both in key order and are checked to hold the same keys. Clusters labelled root (taxid 1: members from different domains) are kept,
so that a near-best hit to one leaves the gene unlabelled; unlabelled clusters (0: no member has a taxon) are left out. FASTA on stdout.
usage: name_reps.py <lca DB> <reps DB> <reps fasta>"""
import sys, os, bisect
lca, reps, fasta = sys.argv[1:]

parts, starts, pos = [], [], 0
i = 0
while os.path.exists(f"{lca}.{i}"):
    starts.append(pos); pos += os.path.getsize(f"{lca}.{i}")
    parts.append(open(f"{lca}.{i}", "rb", buffering=1 << 20)); i += 1
if not parts:
    parts, starts = [open(lca, "rb", buffering=1 << 20)], [0]

def label(off, ln):
    p = bisect.bisect_right(starts, off) - 1
    f = parts[p]
    local = off - starts[p]
    if f.tell() != local:
        f.seek(local)
    return f.read(ln).split(b"\t", 1)[0]

count = {}
out = sys.stdout.buffer
keep = False
with open(lca + ".index", "rb", buffering=1 << 24) as lca_index, open(reps + ".index", "rb", buffering=1 << 24) as rep_index, \
        open(fasta, "rb", buffering=1 << 24) as seqs:
    for line in seqs:
        if line[:1] == b">":
            rep_key = rep_index.readline().split(b"\t", 1)[0]
            k, off, ln = lca_index.readline().split()
            if k != rep_key:
                sys.exit(f"name_reps: key mismatch: representatives {rep_key!r} vs labels {k!r}")
            t = label(int(off), int(ln))
            keep = t != b"0"
            if keep:
                n = count.get(t, 0) + 1; count[t] = n
                out.write(b">" + t + b"." + str(n).encode() + b"\n")
        elif keep:
            out.write(line)
    if rep_index.readline() or lca_index.readline():
        sys.exit("name_reps: index entries left over after the last sequence")
