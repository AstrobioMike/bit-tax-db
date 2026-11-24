#!/usr/bin/env python3
"""Lineage tables for the release: an ID, then names from domain down joined by ";" (the dense taxdump holds only the 7 ranks under root;
skipped ranks are filled, see lineage()).
  write_lineages.py proteins <taxdump dir> <lca tsv>                  > protein-lineages.tsv
      one row per taxon used as a cluster label (column 2 of the LCA table); the ID is the dense taxid; root (1) has an empty lineage
  write_lineages.py genomes <taxdump dir> <taxid.map> <genomes.tsv>   > genome-lineages.tsv
      one row per genome of the k-mer table (column 1 of its genomes.tsv), in that order; the ID is the genome accession"""
import sys
mode, taxdump = sys.argv[1], sys.argv[2]
RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
parent, rank, name = {}, {}, {}
for line in open(f"{taxdump}/nodes.dmp"):
    f = line.split("\t|\t"); parent[f[0]] = f[1]; rank[f[0]] = f[2]
for line in open(f"{taxdump}/names.dmp"):
    f = line.split("\t|\t")
    if f[3].startswith("scientific name"): name[f[0]] = f[1]

def lineage(t):
    """ names from domain down to t's rank; a rank the tree skips (some eukaryotes and viruses) is filled as
    "unclassified <the nearest named taxon above> <rank>", as taxonkit reformat does, so every name sits at its rank's position """
    if t == "1":
        return ""
    named = {}
    while t != "1":
        named[rank[t]] = name[t]; t = parent[t]
    deepest = max(RANKS.index(r) for r in named)
    path, last = [], None
    for r in RANKS[:deepest + 1]:
        if r in named:
            last = named[r]; path.append(last)
        else:
            path.append(f"unclassified {last} {r}")
    return ";".join(path)

out = sys.stdout
if mode == "proteins":
    used = set()
    for line in open(sys.argv[3], buffering=1 << 24):
        used.add(line.split("\t", 2)[1])
    used.discard("0")
    for t in sorted(used, key=int):
        out.write(t + "\t" + lineage(t) + "\n")
elif mode == "genomes":
    taxid = dict(line.rstrip("\n").split("\t")[:2] for line in open(sys.argv[3]))
    for line in open(sys.argv[4]):
        genome = line.split("\t", 1)[0]
        out.write(genome + "\t" + lineage(taxid[genome]) + "\n")
else:
    sys.exit(f"write_lineages: unknown mode {mode!r}")
