#!/usr/bin/env python3
"""Renumber a taxonkit-made taxdump to consecutive taxids (root stays 1).
MMseqs2 sizes its taxonomy tables by the largest taxid, and taxonkit's hashed ids run to 2^31 (an 8.6 GB file)."""
import sys
src, dst = sys.argv[1], sys.argv[2]
ids = sorted(int(l.split("\t|\t")[0]) for l in open(f"{src}/nodes.dmp"))
assert ids[0] == 1
new = {old: i + 1 for i, old in enumerate(ids)}
with open(f"{dst}/nodes.dmp", "w") as out:
    for l in open(f"{src}/nodes.dmp"):
        f = l.split("\t|\t")
        f[0], f[1] = str(new[int(f[0])]), str(new[int(f[1])])
        out.write("\t|\t".join(f))
with open(f"{dst}/names.dmp", "w") as out:
    for l in open(f"{src}/names.dmp"):
        f = l.split("\t|\t")
        f[0] = str(new[int(f[0])])
        out.write("\t|\t".join(f))
with open(f"{dst}/taxid.map", "w") as out:
    for l in open(f"{src}/taxid.map"):
        k, t = l.rstrip("\n").split("\t")
        out.write(f"{k}\t{new[int(t)]}\n")
open(f"{dst}/merged.dmp", "w").close(); open(f"{dst}/delnodes.dmp", "w").close()
print(len(new), "nodes, largest taxid", len(new))
