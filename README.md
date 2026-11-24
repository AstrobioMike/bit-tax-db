# bit-tax-db

This repo holds the reference database that [bit](https://github.com/AstrobioMike/bit) uses for taxonomic classification of contigs and genes (via `classify-seqs` and the `bit run-wf metagenomics` workflow), and the code that builds it.

Current version: **2026-11**, built from GTDB release 232 and NCBI data retrieved 17-Nov-2026.

The database itself is hosted in parts as assets in this repository's [releases](../../releases). Downloading and install is handled by [bit](https://github.com/AstrobioMike/bit) with the programs noted above.


## Database files
Details on making it are below.

### For version 2026-11

| File | What it holds |
|---|---|
| `proteins.dmnd` (117 GB) | DIAMOND database of 353,638,308 protein-cluster representatives, each named `<lineage ID>.<n>` |
| `protein-lineages.tsv` | lineage ID, then the names from domain down to the cluster label's rank, joined by `;` (241,745 IDs) |
| `genome-kmers.hashes.u64` (32 GB) | sampled DNA k-mers of 203,252 genomes (sourmash FracMinHash, k = 31, scaled = 200), one sorted block of uint64 per genome |
| `genome-kmers.genomes.tsv` | genome accession, offset of its block, number of hashes |
| `genome-lineages.tsv` | genome accession, then its names from domain to species |
| `DB-INFO.txt` | sources and versions |

**Sources** 659,705,758 proteins:
- **GTDB R232:** the 199,923 species-representative genomes (618.6 M proteins)
- **NCBI eukaryotes:** 2,811 fungi and 517 protists (one genome per species), plus human GRCh38 (40.3 M proteins)
  - The list is in `data/2026-11/eukaryote-genomes.tsv`; how it was chosen is under [Choosing the eukaryotic genomes](#choosing-the-eukaryotic-genomes).
- **NCBI RefSeq viral proteins:** 15,265 viral taxa (0.7 M proteins).

**Taxonomy.** Bacteria and archaea carry GTDB names; eukaryotes and viruses carry NCBI names. Every lineage has the seven standard ranks merged into one tree with `taxonkit create-taxdump` and renumbered to consecutive IDs.
- **Skipped ranks.** Where a lineage skips a rank (some eukaryotes and many viruses), the rank is filled as
  `unclassified <nearest named taxon above> <rank>`, so every name sits at its rank's position. 
- **Viral entries left out.** Viral proteins are keyed by NCBI taxid. Those whose taxid isn't under Viruses are dropped (7 in 2026-11: 6 keyed
  by their bacterial or archaeal host, and one a cloning vector). 

**Protein clusters.**
- **Clustering.** All proteins were clustered together, GTDB and NCBI in one pass: MMseqs2 `linclust`, 80% identity over 80% of the member.
- **Labels.** Each cluster is labelled with the lowest common ancestor (LCA) of its members' taxa: 78.4% were at species, 17.5% genus, 2.3% family, and 1.8% were higher.
- **Root clusters.** 181,777 clusters (0.05%) mix domains and are labelled root (lineage ID 1, an empty lineage). 
- **Unlabelled clusters.** Clusters with no member taxon are left out of the database.

## How classify-seqs uses it

1. **Species from DNA k-mers.** Each sequence is sketched the same way as the genomes (k = 31, 1 in 200 k-mers). An input sequence gets a species call from this when the best genome shares at least 65% of the input sequence's sampled k-mers (with a minimum requirement of 3 being shared), AND every genome within 90% of the best count agrees at species.
2. **Otherwise, from protein hits.**
   - Genes of an input sequence are called with prodigal (`-p meta`) and searched with DIAMOND blastp (e ≤ 1e-5).
   - **Gene label:** the LCA of the hits within 10% of the best bitscore, capped by the highest identity among those hits: species needs ≥ 95%, genus ≥ 80%, family ≥ 60%; below 60% the label goes no deeper than order.
   - **Input sequence taxonomic assignment:** the deepest taxon holding more than half of its genes' summed bitscores.
   - A k-mer species call takes precedence over the protein label.


## Getting the database

bit downloads and checks it itself (`bit classify-seqs --setup-dbs-only`, or `bit run-wf metagenomics --setup-dbs-only`).

But here's how to do it by hand if wanted: 

For v2026-11, it's 49 parts, 97 GB to download, 157 GB unpacked (about 260 GB free needed while unpacking):

```bash
tag=2026-11; name=bit-tax-db-$tag
base=https://github.com/AstrobioMike/bit-tax-db/releases/download/$tag
curl -sSLO $base/$name.parts.md5
cut -d' ' -f3 $name.parts.md5 | sed "s|^|$base/|" > urls.txt
aria2c -c -x 4 -j 4 -i urls.txt
md5sum -c $name.parts.md5
cat $name.tar.zst.part-* | zstd -dc | tar -xf -
rm $name.tar.zst.part-*
```

## Building a new version (log for future me)

Every version is built the same way, 2026-11 included:

1. **Set the version** at the top of `build.sh`: `VERSION` (year and month, e.g. `2026-11`) and `GTDB_RELEASE` (the GTDB release version wanted as the foundation, e.g., `232.0`.

2. **Build:**

   ```bash
   conda env create -f environment.yml
   conda activate bit-tax-db
   MEM_GB=214 bash build.sh work                # MEM_GB: the memory the job may use
   ```

3. **Write `data/<VERSION>/DB-INFO.txt`** (sources, dates, counts; the counts are in `work/progress.log`). `build.sh` stops before its last
   step (packing) until this file exists. Rerun it to finish; the finished steps are skipped.

4. **Commit** `data/<VERSION>/` and any code changes, and push, so the release's source archive holds the code and lists that built it.

5. **Release** from `work/parts/`, with the GitHub CLI (the tag is the version):

   ```bash
   v=2026-11
   gh release create $v --repo AstrobioMike/bit-tax-db --title "bit-tax-db $v" --notes-file ../release/DB-INFO.txt
   gh release upload $v --repo AstrobioMike/bit-tax-db bit-tax-db-$v.parts.md5 bit-tax-db-$v.tar.zst.part-*
   ```

6. **Point bit at it:** `TAX_DB_VERSION` in bit's `rules/classify-seqs.smk`.

### What build.sh does

It runs every step in order and logs START / DONE / FAIL lines to `work/progress.log`. A rerun skips finished steps, and linclust
resumes from its scratch. It stops a step when free disk falls under `GUARD_GB` (default 30).

**Resources.** 2026-11 was built on 64 cores with 245 GB of RAM and a 2 TB volume:
- **Memory.** Clustering needs the compressed sequences (127 GB) to stay in the page cache next to linclust's own memory (about 41 GB), so plan on at least 256 GB of RAM for the foreseeable future. `build.sh` refuses to start the clustering when `MEM_GB` is too small for that.
- **Disk.** Peak use is during clustering: the compressed sequences, their headers and lookup (79 GB), and linclust's scratch (more than 280 GB).

| Step | What it does | Time for 2026-11 |
|---|---|---|
| 1. gtdb-download | GTDB taxonomy and protein archive (132 GB), md5s from GTDB's `MD5SUM.txt` | 10 min |
| 2. ncbi-download | NCBI taxonomy; the eukaryote list if `data/<VERSION>/` has none; eukaryote proteins (`datasets`); RefSeq viral proteins | minutes |
| 3. ncbi-proteins | NCBI proteins as one FASTA, and `keys.tsv` | minutes |
| 4. taxonomy | GTDB + NCBI lineages → one 7-rank taxonomy (`make_lineages.py`, `fix_lineages.py`, taxonkit, `dense_taxdump.py`) | 1 min |
| 5. createdb | all proteins into one MMseqs2 database, with taxonomy | 1.5 h |
| 6. compress | compressed copy of the sequences | 25 min |
| 7. cluster | linclust 80% / 80%, then each cluster's LCA | ~7 h + 20 min |
| 8. proteins | `protein-lineages.tsv`; representatives renamed `<lineage ID>.<n>` → `proteins.dmnd` | 2.5 h |
| 9. mixed | count of clusters mixing GTDB, eukaryote and viral proteins (a check) | 1 h |
| 10. sketch-gtdb, sketch-ncbi | DNA sketches of all genomes, streamed (nothing kept on disk) | GTDB 1.5 h, NCBI 6.5 h (download-bound) |
| 11. kmers | flat k-mer table and `genome-lineages.tsv` | 3 min |
| 12. package | `release/` → `parts/<name>.tar.zst.part-NNN` + `parts/<name>.parts.md5` | 35 min |

**Scripts** (in `scripts/`, all called by `build.sh`):

| Script | Step | What it does |
|---|---|---|
| `select_eukaryotes.py` | 2 | picks the eukaryote genomes |
| `ncbi_prep.py` | 3 | builds the NCBI FASTA and `keys.tsv` |
| `make_lineages.py` | 4 | writes the merged lineage table |
| `fix_lineages.py` | 4 | guard: drops rows that would give a taxon a second parent |
| `dense_taxdump.py` | 4 | renumbers the taxonomy to consecutive IDs |
| `gtdb_stream.py` | 5 | GTDB archive → FASTA |
| `write_mapping.py` | 5 | protein → taxid mapping |
| `name_reps.py` | 8 | renames the representatives to `<lineage ID>.<n>` |
| `write_lineages.py` | 8, 11 | writes the lineage tables |
| `sketch_gtdb_genomes.py`, `sketch_ncbi_genomes.py`, `sketch_local.py`, `http_reader.py` | 10 | genome sketches |
| `zip2blob.py` | 11 | builds the k-mer table |
| `package.sh` | 12 | packs `release/` into parts |

### Choosing the eukaryotic genomes

One genome with NCBI-annotated proteins per fungal or protist species (protist: a eukaryote that isn't an animal, land plant or fungus),
picked in this order: NCBI-designated reference genome first, then RefSeq over GenBank, then the most protein-coding genes. Human GRCh38 is
added. There is no completeness or assembly-level filter; in 2026-11, 2,364 of the 3,328 picks are reference genomes, 811 are RefSeq, and
1,186 are contig-level assemblies.

## Citing

If you use this databas please cite the sources and tools the database is built from:
- GTDB, NCBI Datasets, NCBI Taxonomy and RefSeq
- MMseqs2, DIAMOND, sourmash, and TaxonKit
