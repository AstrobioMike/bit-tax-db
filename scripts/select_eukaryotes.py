#!/usr/bin/env python3
"""Picks the eukaryote genomes for a new database version (build.sh step 2, when data/<version>/eukaryote-genomes.tsv isn't there yet):
one genome with annotated proteins per fungal or protist species (protist = a eukaryote that isn't an animal, land plant or fungus):
the NCBI reference genome first, then RefSeq over GenBank, then the most protein-coding genes. Human (GRCh38) is added.
Inputs (made in build.sh step 2):
  annotated.tsv   datasets summary genome taxon 2759 --annotated | dataformat tsv genome --fields accession,organism-tax-id,organism-name,
                  assminfo-level,assminfo-refseq-category,source_database,annotinfo-featcount-gene-protein-coding,assmstats-total-sequence-len,
                  assminfo-paired-assm-accession,annotinfo-provider
  taxid-group.tsv taxid, group (fungi / animal / land plant / protist), species name
usage: select_eukaryotes.py annotated.tsv taxid-group.tsv > eukaryote-genomes.tsv"""
import sys, csv, collections
grp, sp = {}, {}
for l in open(sys.argv[2]):
    t, g, s = l.rstrip('\n').split('\t'); grp[t] = g; sp[t] = s
best = {}
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    t = r['Organism Taxonomic ID']; g = grp.get(t)
    if g not in ('fungi', 'protist') or not sp.get(t): continue
    n = int(r['Annotation Count Gene Protein-coding'] or 0)
    if n == 0: continue
    rank = (r['Assembly Refseq Category'] != 'reference genome', not r['Assembly Accession'].startswith('GCF'), -n)
    k = (g, sp[t])
    if k not in best or rank < best[k][0]: best[k] = (rank, r)
o = sys.stdout
o.write("accession\ttaxid\tgroup\tspecies\tproteins\tlevel\trefseq_category\n")
for (g, s), (rk, r) in sorted(best.items()):
    o.write("\t".join([r['Assembly Accession'], r['Organism Taxonomic ID'], g, s, r['Annotation Count Gene Protein-coding'], r['Assembly Level'],
                       r['Assembly Refseq Category']]) + "\n")
o.write("GCF_000001405.40\t9606\thuman\tHomo sapiens\t136807\tChromosome\treference genome\n")
print(collections.Counter(g for g, s in best), file=sys.stderr)
