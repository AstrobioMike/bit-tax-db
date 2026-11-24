#!/usr/bin/env python3
"""The merged lineage table taxonkit create-taxdump builds the taxonomy from: key \t domain..species.
GTDB rows: every genome of the GTDB taxonomy files (ar53 then bac120; RS_/GB_ and the d__ .. s__ prefixes removed).
NCBI rows: ncbi-lineage.tsv (key, NCBI taxid, full NCBI lineage, phylum..species), made by
  cut -f1,2 ncbi/keys.tsv | taxonkit lineage --data-dir ncbi/taxdump -i 2 \\
    | taxonkit reformat2 --data-dir ncbi/taxdump -I 2 -f "{phylum}\\t{class}\\t{order}\\t{family}\\t{genus}\\t{species}"
  domain Viruses for keys starting "v" (viral proteins, keyed by taxid), Eukaryota for the rest (assembly accessions).
  A viral key whose NCBI lineage isn't under Viruses is left out (Nov 22: 6 keyed by their bacterial or archaeal host's taxid, which
  would put whole phyla under Viruses, and 1 artificial cloning vector); these go to stderr.
usage: make_lineages.py ar53_taxonomy.tsv.gz bac120_taxonomy.tsv.gz ncbi-lineage.tsv > lineages.tsv"""
import sys, gzip
ar53, bac120, ncbi = sys.argv[1:]
out = sys.stdout
out.write("key\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
for path in (ar53, bac120):
    with gzip.open(path, "rt") as f:
        for line in f:
            acc, lineage = line.rstrip("\n").split("\t")
            out.write(acc[3:] + "\t" + "\t".join(name[3:] for name in lineage.split(";")) + "\n")
for line in open(ncbi):
    f = line.rstrip("\n").split("\t")
    viral = f[0].startswith("v")
    if viral and not f[2].startswith("Viruses"):
        print(f"left out: {f[0]}\t{f[2]}", file=sys.stderr)
        continue
    out.write(f[0] + "\t" + ("Viruses" if viral else "Eukaryota") + "\t" + "\t".join(f[3:]) + "\n")
