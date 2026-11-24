#!/usr/bin/env python3
"""Reads the GTDB protein_faa_reps tar (uncompressed tar on stdin) and writes one FASTA to stdout,
headers renamed to >{assembly accession}__{original protein id}. Per-genome protein and residue
counts go to the file named as argv[1]. Optional argv[2]: stop after that many genomes (testing)."""
import sys, tarfile, zlib
counts = open(sys.argv[1], "w")
limit = int(sys.argv[2]) if len(sys.argv) > 2 else None
out = sys.stdout.buffer
n = 0
with tarfile.open(fileobj=sys.stdin.buffer, mode="r|") as tar:
    for m in tar:
        if not m.isfile() or not m.name.endswith("_protein.faa.gz"):
            continue
        acc = m.name.rsplit("/", 1)[1][3:-len("_protein.faa.gz")]      # RS_GCF_000016525.1_protein.faa.gz -> GCF_000016525.1
        domain = m.name.split("/")[1]
        data = zlib.decompress(tar.extractfile(m).read(), 31)
        prefix = b">" + acc.encode() + b"__"
        nprot = nres = 0
        lines = []
        for line in data.split(b"\n"):
            if line.startswith(b">"):
                lines.append(prefix + line[1:].split(b" ", 1)[0])
                nprot += 1
            elif line:
                line = line.rstrip(b"*")
                lines.append(line)
                nres += len(line)
        lines.append(b"")
        out.write(b"\n".join(lines))
        counts.write(f"{acc}\t{domain}\t{nprot}\t{nres}\n")
        n += 1
        if limit and n >= limit:
            break
counts.close()
