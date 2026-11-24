#!/usr/bin/env python3
"""Write an MMseqs2 <db>_mapping (internal id, taxid) from <db>.lookup; the taxid comes from the key before '__' in each name.
Keys without a taxon (rows fix_lineages.py dropped) are left out of the mapping, so their proteins count as unclassified in the LCAs."""
import sys, collections
taxmap, lookup, out = sys.argv[1:4]
t = dict(l.rstrip("\n").split("\t") for l in open(taxmap))
n, missing = 0, collections.Counter()
with open(lookup) as fin, open(out, "w", buffering=1 << 24) as fout:
    for l in fin:
        i, name, _ = l.split("\t", 2)
        key = name.partition("__")[0]
        if key in t:
            fout.write(i + "\t" + t[key] + "\n"); n += 1
        else:
            missing[key] += 1
print(n, "entries;", sum(missing.values()), "proteins of", len(missing), "keys without a taxon:", dict(missing))
