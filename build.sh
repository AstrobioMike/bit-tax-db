#!/usr/bin/env bash
# Builds the bit-tax-db reference database from its sources and packs it for a release. See README.md.
#   usage: MEM_GB=<memory for the job, GB> bash build.sh <work dir>
# Run inside the bit-tax-db conda environment (environment.yml). Every step logs START / DONE / FAIL to <work dir>/progress.log and
# leaves a marker in <work dir>/done/; a rerun skips finished steps (and linclust resumes from its scratch).
# Environment: MEM_GB (required): the memory the job may use; THREADS (default 56); GUARD_GB (default 30): stop a step when free disk
# in the work dir falls under this.
set -uo pipefail
REPO=$(cd "$(dirname "$0")" && pwd); S=$REPO/scripts
W=${1:?usage: MEM_GB=<GB> bash build.sh <work dir>}
mkdir -p "$W" && cd "$W" || exit 1
THREADS=${THREADS:-56}
: "${MEM_GB:?set MEM_GB to the memory this job may use, in GB}"
# this database version (the year and month of the build; a second build in the same month is 2026-11b), and the GTDB release built from.
# data/<VERSION>/ holds the version's eukaryote genome list (written by step 2 if it isn't there) and its DB-INFO.txt.
VERSION=2026-11; GTDB_RELEASE=232.0
R=${GTDB_RELEASE%%.*}
NAME=bit-tax-db-$VERSION
DATA=$REPO/data/$VERSION
mkdir -p "$DATA"
GTDB=https://data.ace.uq.edu.au/public/gtdb/data/releases/release$R/$GTDB_RELEASE
mkdir -p dl ncbi tax db clu sketch release parts tmp logs done

log() { echo "$(date -u +%FT%TZ) $*" >> progress.log; }
fail() { log "FAIL $*"; exit 1; }
finished() { [ -e "done/$1" ]; }
finish() { touch "done/$1"; log "DONE $1${2:+: $2}"; }
free_gb() { df --output=avail -BG . | tail -1 | tr -dc 0-9; }
# runs a command in the background and stops it if free disk falls under GUARD_GB
guarded() { "$@" & local pid=$!
  while kill -0 $pid 2>/dev/null; do
    [ "$(free_gb)" -lt "${GUARD_GB:-30}" ] && { log "STOPPING $1: under ${GUARD_GB:-30} GB free disk"; kill $pid; }; sleep 10
  done; wait $pid; }

# 1. GTDB: taxonomy, the representative genomes' proteins (132 GB, removed after step 5); each checked against the md5 in the
#    release's MD5SUM.txt
if ! finished gtdb-download; then
  log "START gtdb-download"
  curl -sSfL "$GTDB/MD5SUM.txt" -o dl/MD5SUM.txt || fail "gtdb-download MD5SUM.txt"
  get() { local md5; md5=$(awk -v f="./$1" '$2 == f {print $1}' dl/MD5SUM.txt); [ -n "$md5" ] || { echo "no md5 for $1"; return 1; }
    aria2c -c -x 16 -s 16 -k 64M --file-allocation=none --console-log-level=warn --summary-interval=0 -d dl --checksum=md5=$md5 "$GTDB/$1"; }
  get ar53_taxonomy_r$R.tsv.gz > logs/gtdb-download.log 2>&1 \
    && get bac120_taxonomy_r$R.tsv.gz >> logs/gtdb-download.log 2>&1 \
    && get genomic_files_reps/gtdb_proteins_aa_reps_r$R.tar.gz >> logs/gtdb-download.log 2>&1 \
    || fail "gtdb-download (logs/gtdb-download.log)"
  finish gtdb-download
fi

# 2. NCBI: taxonomy, the eukaryote genomes' proteins, RefSeq viral proteins.
#    The eukaryote genomes are those in data/<VERSION>/eukaryote-genomes.tsv. If that isn't there it's made from NCBI's current holdings
#    (select_eukaryotes.py: one annotated genome per fungal or protist species, plus human) and written there, to be committed with the
#    version; rebuilding a released version uses its list as kept.
if ! finished ncbi-download; then
  log "START ncbi-download"
  { aria2c -c -x 4 -d ncbi --console-log-level=warn --summary-interval=0 https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
      && mkdir -p ncbi/taxdump && tar -xzf ncbi/new_taxdump.tar.gz -C ncbi/taxdump; } > logs/ncbi-download.log 2>&1 || fail "ncbi-download taxonomy"
  if [ ! -s "$DATA/eukaryote-genomes.tsv" ]; then
    { datasets summary genome taxon 2759 --annotated --as-json-lines \
        | dataformat tsv genome --fields accession,organism-tax-id,organism-name,assminfo-level,assminfo-refseq-category,source_database,annotinfo-featcount-gene-protein-coding,assmstats-total-sequence-len,assminfo-paired-assm-accession,annotinfo-provider \
        > ncbi/euk-annotated.tsv \
      && tail -n +2 ncbi/euk-annotated.tsv | cut -f2 | sort -u | taxonkit lineage --data-dir ncbi/taxdump -t \
        | taxonkit reformat2 --data-dir ncbi/taxdump -I 1 -f "{species}" \
        | awk -F'\t' '{g = "protist"; if ($3 ~ /(^|;)4751(;|$)/) g = "fungi"; else if ($3 ~ /(^|;)33208(;|$)/) g = "animal";
                      else if ($3 ~ /(^|;)3193(;|$)/) g = "land plant"; print $1"\t"g"\t"$4}' > ncbi/euk-taxid-group.tsv \
      && python3 "$S/select_eukaryotes.py" ncbi/euk-annotated.tsv ncbi/euk-taxid-group.tsv > ncbi/eukaryote-genomes.tsv \
      && mv ncbi/eukaryote-genomes.tsv "$DATA/eukaryote-genomes.tsv"; } >> logs/ncbi-download.log 2>&1 || fail "ncbi-download eukaryote selection"
    log "eukaryote genomes selected: $DATA/eukaryote-genomes.tsv ($(tail -n +2 "$DATA/eukaryote-genomes.tsv" | wc -l) genomes)"
  fi
  tail -n +2 "$DATA/eukaryote-genomes.tsv" | cut -f1 > ncbi/eukaryote-accessions.txt
  ( mkdir -p ncbi/euk && cd ncbi/euk \
      && datasets download genome accession --inputfile ../eukaryote-accessions.txt --include protein --dehydrated --filename euk.zip \
      && unzip -q -o euk.zip -d pkg && datasets rehydrate --directory pkg --max-workers 16 ) >> logs/ncbi-download.log 2>&1 \
    || fail "ncbi-download eukaryote proteins"
  mkdir -p ncbi/viral
  curl -sSL https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ | grep -o 'viral\.[0-9]*\.protein\.gpff\.gz' | sort -u \
    | sed 's|^|https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/|' > ncbi/viral/urls.txt
  [ -s ncbi/viral/urls.txt ] && aria2c -c -j 4 -x 4 -d ncbi/viral -i ncbi/viral/urls.txt --console-log-level=warn --summary-interval=0 >> logs/ncbi-download.log 2>&1 \
    || fail "ncbi-download viral proteins"
  finish ncbi-download "$(find ncbi/euk/pkg -name protein.faa | wc -l) eukaryote genomes, $(ls ncbi/viral/*.gpff.gz | wc -l) viral files"
fi

# 3. NCBI proteins as one FASTA (headers <key>__<protein>) + keys.tsv (key, taxid, group, proteins, residues)
if ! finished ncbi-proteins; then
  python3 "$S/ncbi_prep.py" "$DATA/eukaryote-genomes.tsv" ncbi/euk/pkg ncbi/viral ncbi/ncbi.faa ncbi/keys.tsv 2> logs/ncbi-prep.log \
    || fail "ncbi-proteins (logs/ncbi-prep.log)"
  finish ncbi-proteins "$(tail -1 logs/ncbi-prep.log)"
fi

# 4. one taxonomy: GTDB names for bacteria and archaea, NCBI names for eukaryotes and viruses, 7 ranks, consecutive taxids.
#    tax/taxdump-dense/taxid.map maps every key (genome accession, or v<NCBI taxid>) to its taxid.
if ! finished taxonomy; then
  log "START taxonomy"
  cut -f1,2 ncbi/keys.tsv | taxonkit lineage --data-dir ncbi/taxdump -i 2 \
    | taxonkit reformat2 --data-dir ncbi/taxdump -I 2 -f "{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}" > tax/ncbi-lineage.tsv \
    && python3 "$S/make_lineages.py" dl/ar53_taxonomy_r$R.tsv.gz dl/bac120_taxonomy_r$R.tsv.gz tax/ncbi-lineage.tsv > tax/lineages.tsv 2> logs/make-lineages.log \
    && python3 "$S/fix_lineages.py" tax/lineages.tsv > tax/lineages.noheader.tsv 2> logs/fix-lineages.log || fail "taxonomy lineages"
  rm -rf tax/taxdump tax/taxdump-dense && mkdir -p tax/taxdump-dense
  { taxonkit create-taxdump -A 1 -R domain,phylum,class,order,family,genus,species tax/lineages.noheader.tsv -O tax/taxdump --force \
      && python3 "$S/dense_taxdump.py" tax/taxdump tax/taxdump-dense; } > logs/taxdump.log 2>&1 || fail "taxonomy taxdump"
  finish taxonomy "$(grep -c 'left out' logs/make-lineages.log) viral keys left out; $(tail -1 logs/fix-lineages.log)"
fi
D=tax/taxdump-dense

# 5. all proteins in one MMseqs2 database (GTDB streamed from its archive, then NCBI), with the taxonomy attached
if ! finished createdb; then
  log "START createdb"
  rm -f db/all db/all.* db/all_*
  ( set -o pipefail; { pigz -dc dl/gtdb_proteins_aa_reps_r$R.tar.gz | python3 "$S/gtdb_stream.py" db/gtdb-genome-counts.tsv && cat ncbi/ncbi.faa; } \
      | mmseqs createdb stdin db/all --dbtype 1 -v 2 ) > logs/createdb.log 2>&1 || fail "createdb"
  mmseqs createbintaxonomy $D/names.dmp $D/nodes.dmp $D/merged.dmp db/all_taxonomy > logs/bintaxonomy.log 2>&1 \
    && python3 "$S/write_mapping.py" $D/taxid.map db/all.lookup db/all_mapping > logs/mapping.log 2>&1 || fail "createdb taxonomy"
  finish createdb "$(wc -l < db/all.index) proteins ($(wc -l < db/gtdb-genome-counts.tsv) GTDB genomes)"
  rm -f dl/gtdb_proteins_aa_reps_r$R.tar.gz ncbi/ncbi.faa
fi

# 6. compressed copy of the sequences: uncompressed (~214 GB) they plus linclust's own ~41 GB don't fit in 245 GB of memory and
#    linclust's rescoring re-reads them from disk; compressed (~136 GB) they fit. Same IDs, headers and clusters.
if ! finished compress; then
  log "START compress"
  rm -f db/allz db/allz.*
  mmseqs compress db/all db/allz --threads $THREADS -v 2 > logs/compress.log 2>&1 || fail "compress"
  for e in _h _h.index _h.dbtype .lookup .source _mapping _taxonomy; do ln -sf all$e db/allz$e; done
  finish compress "$(du -cb db/allz db/allz.[0-9]* 2>/dev/null | tail -1 | awk '{printf "%.0f GB", $1/1e9}')"
  rm -f db/all db/all.index db/all.dbtype
fi

# 7. clusters at 80% identity over 80% of the member (linclust), then each cluster's label = the LCA of its members' taxa
if ! finished cluster; then
  seq_gb=$(du -cb db/allz db/allz.[0-9]* 2>/dev/null | tail -1 | awk '{printf "%d", $1/1e9 + 0.5}')
  need_gb=$(( seq_gb + 41 + 25 ))   # sequences in the page cache + linclust's own memory + margin
  [ $need_gb -le $MEM_GB ] || fail "cluster not started: needs ~$need_gb GB of memory, MEM_GB is $MEM_GB"
  log "START cluster"
  guarded mmseqs linclust db/allz clu/all80 tmp/linclust --min-seq-id 0.80 -c 0.8 --cov-mode 1 \
      --threads $THREADS --split-memory-limit 100G --remove-tmp-files 1 > logs/linclust.log 2>&1
  [ -s clu/all80.index ] || fail "cluster linclust (its scratch stays in tmp/linclust; a rerun resumes)"
  { mmseqs lca db/allz clu/all80 clu/all80lca --blacklist "" --threads 16 -v 2 \
      && mmseqs createtsv db/allz clu/all80lca clu/all80lca.tsv --threads 16 -v 2; } > logs/lca.log 2>&1 || fail "cluster lca"
  cut -f3 clu/all80lca.tsv | sort -S 4G | uniq -c | sort -k1,1nr > clu/all80lca.ranks.txt
  finish cluster "$(wc -l < clu/all80.index) clusters"
fi

# 8. the protein side of the release: lineage table, and the representatives named <lineage ID>.<n> as a DIAMOND database
if ! finished proteins; then
  log "START proteins"
  python3 "$S/write_lineages.py" proteins $D clu/all80lca.tsv > release/protein-lineages.tsv || fail "proteins lineages"
  rm -f db/rep80 db/rep80.*
  mmseqs createsubdb clu/all80 db/allz db/rep80 --subdb-mode 1 -v 1 > logs/rep80.log 2>&1 \
    && guarded mmseqs convert2fasta db/rep80 tmp/rep80.faa -v 1 >> logs/rep80.log 2>&1 || fail "proteins fasta"
  ( set -o pipefail; python3 "$S/name_reps.py" clu/all80lca db/rep80 tmp/rep80.faa \
      | diamond makedb --db release/proteins.part --threads 16 ) > logs/makedb.log 2>&1 || fail "proteins makedb"
  mv release/proteins.part.dmnd release/proteins.dmnd
  rm -f tmp/rep80.faa
  finish proteins "$(grep -E 'Database sequences' logs/makedb.log | tr -s ' ')"
fi

# 9. how many clusters mix GTDB, NCBI eukaryote and NCBI viral proteins (a check, not part of the release)
if ! finished mixed; then
  mmseqs createtsv db/allz db/allz clu/all80 tmp/all80-members.tsv --threads 8 -v 1 > logs/mixed.log 2>&1 || fail "mixed members"
  awk -F'\t' 'NR == FNR { if ($1 !~ /^v/) euk[$1] = 1; next }
    { split($2, k, "__"); s = (k[1] in euk) ? "E" : ($2 ~ /^GC[AF]_/ ? "G" : "V")
      if ($1 != rep) { if (rep != "") tally(); rep = $1; delete has } has[s] = 1 }
    function tally() { n++; if (("G" in has) && ("E" in has)) ge++; if (("G" in has) && ("V" in has)) gv++; if (("E" in has) && ("V" in has)) ev++ }
    END { tally(); printf "%d clusters; with GTDB + NCBI eukaryote proteins %d, GTDB + NCBI viral %d, eukaryote + viral %d\n", n, ge, gv, ev }' \
    ncbi/keys.tsv tmp/all80-members.tsv > clu/all80-mixed.txt || fail "mixed count"
  rm -f tmp/all80-members.tsv
  finish mixed "$(cat clu/all80-mixed.txt)"
fi

# 10. DNA k-mers of every genome (sourmash FracMinHash, k = 31, scaled = 200), streamed; no genome is kept on disk. NCBI genomes whose
#     stream keeps breaking (large ones, e.g. human) are downloaded with resume, sketched from the file, and removed.
if ! finished sketch-gtdb; then
  log "START sketch-gtdb"
  python3 "$S/sketch_gtdb_genomes.py" $GTDB/genomic_files_reps/gtdb_genomes_reps_r$R.tar.gz sketch/gtdb.sig.zip sketch/gtdb.genome-bp.tsv 16 \
    > logs/sketch-gtdb.log 2>&1 || fail "sketch-gtdb"
  finish sketch-gtdb "$(wc -l < sketch/gtdb.genome-bp.tsv) genomes"
fi
if ! finished sketch-ncbi; then
  log "START sketch-ncbi"
  python3 "$S/sketch_ncbi_genomes.py" ncbi/eukaryote-accessions.txt sketch/ncbi.sig.zip sketch/ncbi.genome-bp.tsv 12 > logs/sketch-ncbi.log 2>&1
  failed=$(grep '^FAILED ' logs/sketch-ncbi.log | cut -d' ' -f2)
  rm -f sketch/ncbi-retry.sig.zip sketch/ncbi-retry.genome-bp.tsv
  if [ -n "$failed" ]; then
    mkdir -p tmp/ncbi-dl
    for acc in $failed; do
      pre=${acc%%_*}; num=${acc#*_}; num=${num%%.*}
      d=https://ftp.ncbi.nlm.nih.gov/genomes/all/$pre/${num:0:3}/${num:3:3}/${num:6:3}/
      dir=$(curl -sSL "$d" | grep -o "href=\"${acc}_[^\"/]*/\"" | head -1 | cut -d'"' -f2 | tr -d /)
      aria2c -c -x 4 -d tmp/ncbi-dl -o $acc.fna.gz --console-log-level=warn --summary-interval=0 "$d$dir/${dir}_genomic.fna.gz" >> logs/sketch-ncbi.log 2>&1 \
        && gzip -t tmp/ncbi-dl/$acc.fna.gz || fail "sketch-ncbi download $acc"
    done
    python3 "$S/sketch_local.py" sketch/ncbi-retry.sig.zip sketch/ncbi-retry.genome-bp.tsv tmp/ncbi-dl/*.fna.gz >> logs/sketch-ncbi.log 2>&1 \
      || fail "sketch-ncbi from files"
    rm -rf tmp/ncbi-dl
  fi
  n=$(( $(wc -l < sketch/ncbi.genome-bp.tsv) + $(cat sketch/ncbi-retry.genome-bp.tsv 2>/dev/null | wc -l) ))
  [ $n -eq $(wc -l < ncbi/eukaryote-accessions.txt) ] || fail "sketch-ncbi: $n of $(wc -l < ncbi/eukaryote-accessions.txt) genomes sketched"
  finish sketch-ncbi "$n genomes ($(echo $failed | wc -w) from downloaded files)"
fi

# 11. the k-mer side of the release: one flat table of all genomes' k-mers, and the genomes' lineages (in the table's order)
if ! finished kmers; then
  log "START kmers"
  zips="sketch/gtdb.sig.zip sketch/ncbi.sig.zip"; [ -s sketch/ncbi-retry.sig.zip ] && zips="$zips sketch/ncbi-retry.sig.zip"
  python3 "$S/zip2blob.py" release/genome-kmers 12 $zips > logs/kmers.log 2>&1 || fail "kmers table"
  python3 "$S/write_lineages.py" genomes $D $D/taxid.map release/genome-kmers.genomes.tsv > release/genome-lineages.tsv || fail "kmers lineages"
  finish kmers "$(tail -1 logs/kmers.log)"
fi

# 12. pack: release/ (+ data/<VERSION>/DB-INFO.txt) -> parts/<name>.tar.zst.part-NNN (each under 2 GiB) + parts/<name>.parts.md5.
#     Stops here until DB-INFO.txt is written (its counts are in progress.log); rerun to pack.
if ! finished package; then
  [ -s "$DATA/DB-INFO.txt" ] || fail "package: write $DATA/DB-INFO.txt (sources, dates, counts from progress.log), then rerun"
  log "START package"
  cp "$DATA/DB-INFO.txt" release/DB-INFO.txt
  bash "$S/package.sh" release parts $NAME > logs/package.log 2>&1 || fail "package"
  finish package "$(tail -1 logs/package.log)"
fi
log "DONE build $NAME"
