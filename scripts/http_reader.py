"""Sequential file-like reader over parallel HTTP range requests, for streaming a large archive without storing it."""
import io, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor


class Http(io.RawIOBase):
    """Reads CONN chunks of CHUNK bytes at a time, in order; each chunk is retried until it arrives whole."""
    CHUNK, CONN = 64 << 20, 16

    def __init__(self, url):
        self.url = url
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "python-urllib"})
        self.size = int(urllib.request.urlopen(req, timeout=120).headers["Content-Length"])
        self.pool = ThreadPoolExecutor(self.CONN)
        self.next, self.pending, self.cur, self.off = 0, [], b"", 0

    def _get(self, start):
        end = min(start + self.CHUNK, self.size) - 1
        for attempt in range(500):
            try:
                req = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{end}", "User-Agent": "python-urllib"})
                d = urllib.request.urlopen(req, timeout=120).read()
                if len(d) == end - start + 1:
                    return d
                raise IOError(f"short read {len(d)}")
            except Exception as e:
                print(f"retry chunk at {start}: {e}", file=sys.stderr, flush=True)
                time.sleep(min(60, 5 * (attempt + 1)))
        raise IOError("gave up")

    def _fill(self):
        while len(self.pending) < self.CONN and self.next < self.size:
            self.pending.append(self.pool.submit(self._get, self.next))
            self.next += self.CHUNK

    def readable(self):
        return True

    def readinto(self, b):
        if self.off >= len(self.cur):
            self._fill()
            if not self.pending:
                return 0
            self.cur, self.off = self.pending.pop(0).result(), 0
        n = min(len(b), len(self.cur) - self.off)
        b[:n] = self.cur[self.off:self.off + n]
        self.off += n
        return n
