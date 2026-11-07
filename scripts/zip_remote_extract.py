"""Extract selected members from a remote ZIP (HTTP range requests; supports ZIP64), verifying each member's
CRC32 against the archive's central directory. Avoids downloading multi-GB bundles.
Usage: python scripts/zip_remote_extract.py <url> <dest_dir> <member_suffix> [<member_suffix> ...]
       python scripts/zip_remote_extract.py <url> --list
"""
import json, os, re, struct, subprocess, sys, time, zlib, hashlib


def rng_get(url, start, end):
    for attempt in range(4):
        r = subprocess.run(["curl", "-sfL", "-r", f"{start}-{end}", url], capture_output=True, timeout=900)
        if r.returncode == 0 and len(r.stdout) == end - start + 1:
            return r.stdout
        print("retry", attempt, r.returncode, len(r.stdout)); time.sleep(5 * (attempt + 1))
    raise RuntimeError("range request failed")


def size_of(url):
    r = subprocess.run(["curl", "-sIL", url], capture_output=True, text=True, timeout=60)
    return int(re.findall(r"(?i)content-length:\s*(\d+)", r.stdout)[-1])


def central_directory(url):
    total = size_of(url)
    tail = rng_get(url, max(0, total - 70000), total - 1)
    i = tail.rfind(b"PK\x05\x06")
    cd_size, cd_off = struct.unpack("<II", tail[i + 12:i + 20])
    n_entries = struct.unpack("<H", tail[i + 10:i + 12])[0]
    if cd_off == 0xFFFFFFFF or cd_size == 0xFFFFFFFF or n_entries == 0xFFFF:   # ZIP64
        j = tail.rfind(b"PK\x06\x07")
        z64_off = struct.unpack("<Q", tail[j + 8:j + 16])[0]
        z = rng_get(url, z64_off, z64_off + 55)
        n_entries, cd_size, cd_off = struct.unpack("<QQQ", z[32:56])
    cd = rng_get(url, cd_off, cd_off + cd_size - 1)
    out, p = [], 0
    while p < len(cd) and cd[p:p + 4] == b"PK\x01\x02":
        (method, crc, csize, usize) = struct.unpack("<HxxxxIII", cd[p + 10:p + 28])
        nlen, elen, clen = struct.unpack("<HHH", cd[p + 28:p + 34])
        loff = struct.unpack("<I", cd[p + 42:p + 46])[0]
        name = cd[p + 46:p + 46 + nlen].decode("utf-8", "replace")
        extra = cd[p + 46 + nlen:p + 46 + nlen + elen]
        q = 0
        while q < len(extra):
            hid, hsz = struct.unpack("<HH", extra[q:q + 4])
            if hid == 1:
                vals = extra[q + 4:q + 4 + hsz]; k = 0
                if usize == 0xFFFFFFFF: usize = struct.unpack("<Q", vals[k:k + 8])[0]; k += 8
                if csize == 0xFFFFFFFF: csize = struct.unpack("<Q", vals[k:k + 8])[0]; k += 8
                if loff == 0xFFFFFFFF: loff = struct.unpack("<Q", vals[k:k + 8])[0]; k += 8
            q += 4 + hsz
        out.append(dict(name=name, method=method, crc=crc, csize=csize, usize=usize, offset=loff))
        p += 46 + nlen + elen + clen
    return total, out


def extract(url, entry, dest):
    h = rng_get(url, entry["offset"], entry["offset"] + 29)
    nlen, elen = struct.unpack("<HH", h[26:30])
    start = entry["offset"] + 30 + nlen + elen
    data = bytearray()
    chunk = 64 << 20
    for s in range(start, start + entry["csize"], chunk):
        data += rng_get(url, s, min(s + chunk, start + entry["csize"]) - 1)
    if entry["method"] == 8:
        raw = zlib.decompress(bytes(data), -15)
    elif entry["method"] == 0:
        raw = bytes(data)
    else:
        raise ValueError(f"unsupported method {entry['method']}")
    assert len(raw) == entry["usize"], "size mismatch"
    assert (zlib.crc32(raw) & 0xFFFFFFFF) == entry["crc"], "CRC32 mismatch"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "wb") as f:
        f.write(raw)
    return hashlib.sha256(raw).hexdigest()


if __name__ == "__main__":
    url = sys.argv[1]
    total, entries = central_directory(url)
    if sys.argv[2] == "--list":
        for e in entries:
            print(e["usize"], e["name"])
        sys.exit()
    dest_dir, wants = sys.argv[2], sys.argv[3:]
    for w in wants:
        cands = [e for e in entries if e["name"].endswith(w)]
        assert len(cands) == 1, (w, [e["name"] for e in cands])
        e = cands[0]
        dest = os.path.join(dest_dir, os.path.basename(e["name"]))
        if os.path.exists(dest) and os.path.getsize(dest) == e["usize"]:
            print("exists", dest); continue
        sha = extract(url, e, dest)
        rec = dict(url=url, zip_bytes=total, member=e["name"], path=dest, bytes=e["usize"],
                   zip_crc32=f"{e['crc']:08x}", crc32_verified=True, sha256=sha,
                   downloaded=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
        with open("data/raw/MANIFEST.jsonl", "a") as m:
            m.write(json.dumps(rec) + "\n")
        print("ok", dest, e["usize"])
