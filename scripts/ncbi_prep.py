#!/usr/bin/env python3
"""NCBI side of the database: one FASTA of the eukaryote genomes' proteins (headers <assembly>__<protein>) and the RefSeq viral
proteins (headers v<taxid>__<protein>), and keys.tsv: key (the text before "__"), NCBI taxid, group, proteins, residues.
usage: ncbi_prep.py eukaryote-genomes.tsv <datasets package dir> <viral gpff dir> out.faa keys.tsv"""
import sys, csv, glob, gzip, os, collections
genomes, pkg, viral, faa, keys_path = sys.argv[1:]
sel = {r["accession"]: r for r in csv.DictReader(open(genomes), delimiter="\t")}
out, keys = open(faa, "w"), open(keys_path, "w")
missing = []
for acc, r in sorted(sel.items()):
    path = f"{pkg}/ncbi_dataset/data/{acc}/protein.faa"
    if not os.path.exists(path):
        missing.append(acc); continue
    nprot = nres = 0
    for line in open(path):
        if line.startswith(">"):
            out.write(f">{acc}__{line[1:].split()[0]}\n"); nprot += 1
        else:
            out.write(line); nres += len(line) - 1
    keys.write(f"{acc}\t{r['taxid']}\t{r['group']}\t{nprot}\t{nres}\n")
vir = collections.defaultdict(lambda: [0, 0])
for path in sorted(glob.glob(f"{viral}/viral.*.protein.gpff.gz")):
    ver = tax = None; seq = []; in_seq = False
    for line in gzip.open(path, "rt"):
        if line.startswith("VERSION"):
            ver = line.split()[1]
        elif '/db_xref="taxon:' in line and tax is None:
            tax = line.split('taxon:')[1].split('"')[0]
        elif line.startswith("ORIGIN"):
            in_seq = True
        elif line.startswith("//"):
            s = "".join(seq).upper()
            if ver and tax and s:
                out.write(f">v{tax}__{ver}\n{s}\n"); vir[tax][0] += 1; vir[tax][1] += len(s)
            ver = tax = None; seq = []; in_seq = False
        elif in_seq:
            seq.append("".join(line.split()[1:]))
for tax, (n, r) in sorted(vir.items()):
    keys.write(f"v{tax}\t{tax}\tvirus\t{n}\t{r}\n")
print(f"{len(sel) - len(missing)} eukaryote genomes, {len(vir)} viral taxa", file=sys.stderr)
if missing:
    sys.exit(f"ncbi_prep: no proteins for {len(missing)} genomes: {' '.join(missing[:20])}")
