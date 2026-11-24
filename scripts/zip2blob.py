#!/usr/bin/env python3
"""Convert sourmash signature zips (one k-mer size, one scaled) into one flat k-mer table for fast scanning:
  <out>.hashes.u64   all genomes' hashes, each genome's block sorted (raw little-endian uint64)
  <out>.genomes.tsv  name, offset (in hashes), count        (in block order)
usage: zip2blob.py out_prefix workers a.sig.zip [b.sig.zip ...]   (duplicate names across zips are taken once)"""
import sys, os, zipfile, gzip, json, numpy as np
from multiprocessing import Pool
out, workers, zips = sys.argv[1], int(sys.argv[2]), sys.argv[3:]

def work(job):
    zi, part, names = job
    z = zipfile.ZipFile(zips[zi]); rows = []
    with open(f"{out}.part{zi}.{part}.u64", "wb") as fh:
        for n in names:
            raw = z.read(n)
            if raw[:2] == b"\x1f\x8b": raw = gzip.decompress(raw)
            for rec in json.loads(raw):
                for s in rec["signatures"]:
                    a = np.array(s["mins"], dtype=np.uint64); a.sort(); a.tofile(fh)
                    rows.append((rec["name"].split()[0], len(a), s["ksize"], 18446744073709551615 // s["max_hash"] if s.get("max_hash") else 0))
    return zi, part, rows

def main():
    jobs = []
    for zi, zp in enumerate(zips):
        names = [n for n in zipfile.ZipFile(zp).namelist() if n.startswith("signatures/")]
        step = max(1, len(names) // (workers * 8) + 1)
        jobs += [(zi, p, names[i:i + step]) for p, i in enumerate(range(0, len(names), step))]
    res = sorted(Pool(workers).imap_unordered(work, jobs, chunksize=1))
    seen, off = set(), 0
    with open(out + ".hashes.u64", "wb") as fh, open(out + ".genomes.tsv", "w") as tsv:
        for zi, part, rows in res:
            pf = f"{out}.part{zi}.{part}.u64"; data = np.fromfile(pf, dtype=np.uint64); p = 0
            for name, n, k, scaled in rows:
                if name not in seen:
                    seen.add(name); data[p:p + n].tofile(fh); tsv.write(f"{name}\t{off}\t{n}\n"); off += n
                p += n
            os.remove(pf)
    print(f"{len(seen)} genomes, {off} hashes, k/scaled of last: {k}/{scaled}")

if __name__ == "__main__":
    main()
