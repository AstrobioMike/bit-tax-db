#!/usr/bin/env python3
"""Stream GTDB's representative-genome archive and write one DNA sketch per genome (k=31, scaled=200), keeping no genome on disk.
The HTTP reader reconnects at its current offset if the connection drops. Output: sourmash zip + a tsv of genome sizes."""
import sys, re, io, gzip, tarfile, time, urllib.request, multiprocessing as mp
import sourmash, screed
from sourmash import MinHash, SourmashSignature
from sourmash.sourmash_args import SaveSignaturesToLocation

URL, OUT, TSV = sys.argv[1], sys.argv[2], sys.argv[3]
K, SCALED, WORKERS = 31, 200, int(sys.argv[4])

from http_reader import Http

def sketch(job):
    acc, data = job
    mh = MinHash(n=0, ksize=K, scaled=SCALED)
    bp = 0
    with gzip.open(io.BytesIO(data), "rt") as f:
        seq = []
        for line in f:
            if line.startswith(">"):
                if seq: s = "".join(seq); bp += len(s); mh.add_sequence(s, force=True); seq = []
            else: seq.append(line.strip())
        if seq: s = "".join(seq); bp += len(s); mh.add_sequence(s, force=True)
    return acc, bp, sourmash.save_signatures_to_json([SourmashSignature(mh, name=acc)])

def members(url):
    tf = tarfile.open(fileobj=io.BufferedReader(Http(url), 1 << 24), mode="r|gz")
    for m in tf:
        if not m.isfile(): continue
        a = re.search(r"(GC[AF]_\d+\.\d+)", m.name)
        if a and m.name.endswith(".gz"):
            yield a.group(1), tf.extractfile(m).read()

if __name__ == "__main__":
    n = 0; t0 = time.time()
    with mp.Pool(WORKERS) as pool, SaveSignaturesToLocation(OUT) as save, open(TSV, "w") as tsv:
        import threading
        slots = threading.BoundedSemaphore(WORKERS * 8)   # backpressure: imap's feeder would otherwise queue the whole archive in memory
        def bounded():
            for job in members(URL):
                slots.acquire(); yield job
        for acc, bp, js in pool.imap_unordered(sketch, bounded(), chunksize=1):
            slots.release()
            for ss in sourmash.load_signatures_from_json(js): save.add(ss)
            tsv.write(f"{acc}\t{bp}\n"); n += 1
            if n % 2000 == 0: print(f"{n} genomes, {time.time()-t0:.0f} s", file=sys.stderr, flush=True); tsv.flush()
    print(f"DONE {n} genomes, {time.time()-t0:.0f} s", file=sys.stderr, flush=True)
