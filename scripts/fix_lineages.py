#!/usr/bin/env python3
"""Drops the rows of the merged lineage table that would give a taxon a second parent, before taxonkit create-taxdump (which keys taxa by
rank + name, so one stray row moves a whole taxon: Nov 22, 6 NCBI viral-protein keys whose taxids are bacterial or archaeal hosts, domain
"Viruses", put all of Bacillota, Methanobacteriota and Mycoplasmatota under Viruses). For each taxon with more than one parent the parent
most rows give is kept; rows naming another parent are dropped (their proteins get no taxon and are left out of the cluster LCAs).
The dropped rows go to stderr. usage: fix_lineages.py lineages.tsv (with header) > lineages.fixed.noheader.tsv"""
import sys, collections
RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
rows = [line.rstrip("\n").split("\t") for line in open(sys.argv[1])][1:]

def edges(row):
    prev = "root"
    for rank, name in zip(RANKS, row[1:]):
        if name:
            yield (rank, name), prev
            prev = name

parents = collections.defaultdict(collections.Counter)
for row in rows:
    for taxon, parent in edges(row):
        parents[taxon][parent] += 1
keep = {taxon: c.most_common(1)[0][0] for taxon, c in parents.items()}
out, dropped = sys.stdout, 0
for row in rows:
    if all(keep[taxon] == parent for taxon, parent in edges(row)):
        out.write("\t".join(row) + "\n")
    else:
        dropped += 1
        print("dropped:", "\t".join(row), file=sys.stderr)
print(f"{dropped} of {len(rows)} rows dropped; {sum(len(c) > 1 for c in parents.values())} taxa had more than one parent", file=sys.stderr)
