#!/usr/bin/env python3
"""DNA sketch (k=31, scaled=200) of local genome files named <accession>.fna.gz. usage: sketch_local.py out.sig.zip sizes.tsv file...
Large genomes: sourmash keeps a sketch as a sorted array, so adding k-mers to a sketch of millions of hashes is quadratic. Each 5 Mb piece is
sketched on its own (pieces overlap by k-1 so no k-mer is lost) and the hashes are combined once at the end."""
import sys, os, gzip
import numpy as np
from sourmash import MinHash, SourmashSignature
from sourmash.sourmash_args import SaveSignaturesToLocation
K, SCALED, PIECE = 31, 200, 5000000
def hashes(s, out):
    for i in range(0, max(1, len(s) - K + 1), PIECE):
        m = MinHash(n=0, ksize=K, scaled=SCALED); m.add_sequence(s[i:i + PIECE + K - 1], force=True)
        out.append(np.fromiter(m.hashes.keys(), dtype=np.uint64))
with SaveSignaturesToLocation(sys.argv[1]) as save, open(sys.argv[2], "a") as tsv:
    for path in sys.argv[3:]:
        acc = os.path.basename(path).replace(".fna.gz", "")
        parts = []; bp = 0; seq = []
        with gzip.open(path, "rt") as f:
            for line in f:
                if line[0] == ">":
                    if seq: s = "".join(seq); bp += len(s); hashes(s, parts); seq = []
                else: seq.append(line.strip())
            if seq: s = "".join(seq); bp += len(s); hashes(s, parts)
        allh = np.unique(np.concatenate(parts))
        mh = MinHash(n=0, ksize=K, scaled=SCALED); mh.add_many(allh.tolist())
        save.add(SourmashSignature(mh, name=acc)); tsv.write(f"{acc}\t{bp}\t{len(allh)}\n"); tsv.flush()
        print(f"{acc}\t{bp}\t{len(allh)}", file=sys.stderr, flush=True)
