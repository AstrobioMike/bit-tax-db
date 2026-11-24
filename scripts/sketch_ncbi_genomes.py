#!/usr/bin/env python3
"""DNA sketch (k=31, scaled=200) of each NCBI assembly in a list, streamed from the NCBI FTP site; no genome is stored.
usage: sketch_ncbi_genomes.py accessions.txt out.sig.zip sizes.tsv workers"""
import sys, re, gzip, time, urllib.request, multiprocessing as mp
import sourmash
from sourmash import MinHash, SourmashSignature
from sourmash.sourmash_args import SaveSignaturesToLocation

K, SCALED = 31, 200
BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

def get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "python-urllib"}), timeout=300)

def sketch_once(acc):
    pre, num = acc.split("_"); num = num.split(".")[0]
    d = f"{BASE}/{pre}/{num[0:3]}/{num[3:6]}/{num[6:9]}/"
    names = re.findall(r'href="(%s_[^"/]+)/"' % re.escape(acc), get(d).read().decode())
    if not names: raise IOError("no directory for " + acc)
    mh = MinHash(n=0, ksize=K, scaled=SCALED); bp = 0; seq = []
    with gzip.open(get(f"{d}{names[0]}/{names[0]}_genomic.fna.gz"), "rt") as f:
        for line in f:
            if line[0] == ">":
                if seq: s = "".join(seq); bp += len(s); mh.add_sequence(s, force=True); seq = []
            else: seq.append(line.strip())
        if seq: s = "".join(seq); bp += len(s); mh.add_sequence(s, force=True)
    return acc, bp, sourmash.save_signatures_to_json([SourmashSignature(mh, name=acc)])

def sketch(acc):
    for attempt in range(8):
        try: return sketch_once(acc)
        except Exception as e:
            print(f"retry {acc}: {e}", file=sys.stderr, flush=True); time.sleep(10 * (attempt + 1))
    return acc, -1, None

if __name__ == "__main__":
    accs = [l.split()[0] for l in open(sys.argv[1]) if l.startswith("GC")]
    n = bad = 0; t0 = time.time()
    with mp.Pool(int(sys.argv[4])) as pool, SaveSignaturesToLocation(sys.argv[2]) as save, open(sys.argv[3], "w") as tsv:
        for acc, bp, js in pool.imap_unordered(sketch, accs, chunksize=1):
            if js is None: bad += 1; print(f"FAILED {acc}", file=sys.stderr, flush=True); continue
            for ss in sourmash.load_signatures_from_json(js): save.add(ss)
            tsv.write(f"{acc}\t{bp}\n"); tsv.flush(); n += 1
            if n % 100 == 0: print(f"{n} genomes, {time.time()-t0:.0f} s", file=sys.stderr, flush=True)
    print(f"DONE {n} genomes, {bad} failed, {time.time()-t0:.0f} s", file=sys.stderr, flush=True)
    sys.exit(1 if bad else 0)
