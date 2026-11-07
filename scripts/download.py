"""Download public files with size verification + sha256 manifest. Never overwrites existing verified files."""
import hashlib, json, os, sys, time, urllib.request
MANIFEST = "data/raw/MANIFEST.jsonl"

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def head_size(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", -1))

def fetch(url, dest, source_page=""):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    expected = head_size(url)
    if os.path.exists(dest) and os.path.getsize(dest) == expected:
        print("exists", dest); return
    tmp = dest + ".part"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b: break
                    f.write(b)
            break
        except Exception as e:
            print("retry", attempt, e); time.sleep(10 * (attempt + 1))
    size = os.path.getsize(tmp)
    assert size == expected, f"size mismatch {size} != {expected} for {url}"
    os.replace(tmp, dest)
    rec = dict(url=url, path=dest, bytes=size, sha256=sha256(dest), source_page=source_page,
               downloaded=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    with open(MANIFEST, "a") as m:
        m.write(json.dumps(rec) + "\n")
    print("ok", dest, size)

if __name__ == "__main__":
    spec = json.load(open(sys.argv[1]))
    for item in spec:
        fetch(item["url"], item["dest"], item.get("source_page", ""))
