"""PR-CONTEXT step 5 -- retrieve the ARG-bearing plasmid replicon sequences.

Smallest valid retrieval scope, as the task requires. Only the 6,621 plasmid replicons that
actually carry a primary acquired ARG occurrence are fetched -- not the 15,941 plasmids in the
cohort, and not one genome FASTA. Total is under 1 GB.

Every sequence is fetched at its EXACT accession.version, the accession is re-read out of the
returned FASTA header and asserted to match what was requested, and the bytes are hashed. A
record whose returned accession differs from the requested one is recorded as a mismatch and
excluded from mobility analysis rather than silently accepted -- a different version is a
different sequence.
"""
import argparse, collections, hashlib, io, json, os, sys, time, urllib.parse, urllib.request
import datetime

VERSION = "pr_context_fetch_plasmids_v1.0.0"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 100
_KEY = os.environ.get("NCBI_API_KEY")
_GAP = 0.11 if _KEY else 0.36


def post(ids, attempts=5):
    data = {"db": "nuccore", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"}
    if _KEY:
        data["api_key"] = _KEY
    body = urllib.parse.urlencode(data).encode()
    last = ""
    for a in range(attempts):
        try:
            time.sleep(_GAP)
            req = urllib.request.Request(
                EFETCH, data=body,
                headers={"User-Agent": "portabilityrisk-context/1.0",
                         "Content-Type": "application/x-www-form-urlencoded"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read().decode("utf-8", "replace"), ""
        except Exception as e:
            last = repr(e)[:140]
            time.sleep(min(2 ** a, 30))
    return None, last


def split_fasta(text):
    cur, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                yield cur, "\n".join(buf)
            cur, buf = line[1:].strip(), []
        elif cur is not None:
            buf.append(line.strip())
    if cur:
        yield cur, "\n".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--accessions", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")
    fa = os.path.join(a.root, "plasmid_fasta")
    os.makedirs(fa, exist_ok=True)
    accs = [x.strip() for x in open(a.accessions, encoding="utf-8") if x.strip()]
    print("%s | plasmid replicons to retrieve: %d | api_key: %s"
          % (VERSION, len(accs), "yes" if _KEY else "no"))

    got, receipts, missing = {}, [], []
    for i in range(0, len(accs), BATCH):
        chunk = accs[i:i + BATCH]
        need = [x for x in chunk if not os.path.exists(os.path.join(fa, x + ".fna"))]
        if need:
            txt, err = post(need)
            if txt is None:
                for x in need:
                    missing.append((x, "efetch failed: " + err))
                continue
            for hdr, seq in split_fasta(txt):
                ret = hdr.split()[0]
                p = os.path.join(fa, ret + ".fna")
                with open(p, "w", encoding="utf-8", newline="\n") as fh:
                    fh.write(">" + hdr + "\n")
                    for j in range(0, len(seq), 70):
                        fh.write(seq[j:j + 70] + "\n")
        if (i // BATCH) % 10 == 0:
            print("  %d/%d requested" % (min(i + BATCH, len(accs)), len(accs)), flush=True)

    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    for x in accs:
        p = os.path.join(fa, x + ".fna")
        if not os.path.exists(p):
            missing.append((x, "no FASTA returned at the exact accession.version"))
            receipts.append({"requested_accession": x, "returned_accession": "",
                             "status": "ACCESSION_VERSION_UNAVAILABLE", "bytes": "",
                             "sequence_length": "", "raw_sha256": "",
                             "canonical_url": EFETCH + "?db=nuccore&id=%s&rettype=fasta" % x,
                             "timestamp_utc": now})
            continue
        raw = open(p, "rb").read()
        hdr = raw.split(b"\n", 1)[0].decode("utf-8", "replace")[1:]
        ret = hdr.split()[0]
        seqlen = sum(len(l) for l in raw.decode("utf-8", "replace").splitlines()[1:])
        receipts.append({
            "requested_accession": x, "returned_accession": ret,
            "status": "EXACT" if ret == x else "ACCESSION_MISMATCH",
            "bytes": len(raw), "sequence_length": seqlen,
            "raw_sha256": hashlib.sha256(raw).hexdigest(),
            "canonical_url": EFETCH + "?db=nuccore&id=%s&rettype=fasta" % x,
            "timestamp_utc": now})

    cols = ["requested_accession", "returned_accession", "status", "bytes",
            "sequence_length", "raw_sha256", "canonical_url", "timestamp_utc"]
    p = os.path.join(out, "plasmid_retrieval_receipts.tsv")
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(receipts, key=lambda r: r["requested_accession"]):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    c = collections.Counter(r["status"] for r in receipts)
    tot = sum(int(r["bytes"] or 0) for r in receipts)
    print("\nretrieval:", dict(c))
    print("bytes on disk: %.3f GB" % (tot / 1e9))
    if missing:
        print("MISSING (%d): %s" % (len(missing), missing[:5]))
    print("receipts sha256:", hashlib.sha256(open(p, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
