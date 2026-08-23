"""PR-CONTEXT step 6 -- extract ARG neighbourhood windows on CHROMOSOMES.

Only the windows are retrieved, never whole chromosomes. Fetching the 6,190 ARG-bearing
chromosomes would be about 31 GB; the merged +/-10 kb windows are about 0.51 GB, and they
contain everything the neighbourhood analysis can use. Plasmid neighbourhoods need no
retrieval at all because the whole plasmid is already held.

Windows that overlap on one replicon are MERGED before retrieval, so a region shared by
several ARGs is downloaded once and recorded once as a shared context block. Each ARG keeps
its own window record pointing into that block, so co-located ARGs are never silently
collapsed into a single observation, and never double-counted as independent contexts.

Circularity is respected: on a circular replicon a window that runs off either end WRAPS and
is recorded as wrapped. Truncation is recorded only where wrapping is not valid.
"""
import argparse, collections, csv, hashlib, json, os, sys, time, urllib.parse, urllib.request
import datetime

VERSION = "pr_context_fetch_windows_v1.0.0"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_KEY = os.environ.get("NCBI_API_KEY")
_GAP = 0.11 if _KEY else 0.36
WIN = 10000


def _post(url, data, attempts=5):
    if _KEY:
        data = dict(data, api_key=_KEY)
    body = urllib.parse.urlencode(data).encode()
    last = ""
    for a in range(attempts):
        try:
            time.sleep(_GAP)
            req = urllib.request.Request(
                url, data=body, headers={"User-Agent": "portabilityrisk-context/1.0"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read().decode("utf-8", "replace"), ""
        except Exception as e:
            last = repr(e)[:140]
            time.sleep(min(2 ** a, 30))
    return None, last


def topology(accs, cache):
    """Circular or linear, from nuccore esummary. Batched 200 at a time. Never assumed."""
    out = {}
    todo = []
    for a in accs:
        p = os.path.join(cache, a + ".topo")
        if os.path.exists(p):
            out[a] = open(p, encoding="utf-8").read().strip()
        else:
            todo.append(a)
    for i in range(0, len(todo), 200):
        chunk = todo[i:i + 200]
        txt, err = _post(ESUMMARY, {"db": "nuccore", "id": ",".join(chunk),
                                    "retmode": "json"})
        if txt is None:
            for a in chunk:
                out[a] = "unknown"
            continue
        try:
            d = json.loads(txt).get("result", {})
        except Exception:
            d = {}
        byacc = {}
        for k, v in d.items():
            if k == "uids" or not isinstance(v, dict):
                continue
            byacc[v.get("accessionversion") or v.get("caption", "")] = (
                v.get("topology") or "unknown")
        for a in chunk:
            t = byacc.get(a, "unknown")
            out[a] = t
            open(os.path.join(cache, a + ".topo"), "w", encoding="utf-8").write(t)
        print("  topology %d/%d" % (min(i + 200, len(todo)), len(todo)), flush=True)
    return out


def fetch_range(acc, start, stop, strand=1):
    txt, err = _post(EFETCH, {"db": "nuccore", "id": acc, "rettype": "fasta",
                              "retmode": "text", "seq_start": str(start),
                              "seq_stop": str(stop), "strand": str(strand)})
    if txt is None:
        return None, err
    lines = txt.splitlines()
    return "".join(l.strip() for l in lines[1:] if not l.startswith(">")), ""


def write(path, rows, cols):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")
    seqdir = os.path.join(a.root, "window_fasta")
    cache = os.path.join(a.root, "topo_cache")
    for d in (seqdir, cache):
        os.makedirs(d, exist_ok=True)

    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"
        and r["evidence_type"] == "direct_chromosome"]
    print("%s | chromosomal ARG occurrences: %d" % (VERSION, len(occ)))
    accs = sorted({r["replicon_accession"] for r in occ})
    print("ARG-bearing chromosomes: %d" % len(accs))
    topo = topology(accs, cache)
    tc = collections.Counter(topo.values())
    print("topology:", dict(tc))

    lens = {r["replicon_accession"]: int(r["replicon_length"]) for r in occ
            if r["replicon_length"]}
    by = collections.defaultdict(list)
    for r in occ:
        by[r["replicon_accession"]].append(r)

    blocks, wrows, failures = [], [], []
    bi = 0
    for k, rows in sorted(by.items()):
        L = lens.get(k, 0)
        circ = topo.get(k, "unknown") == "circular"
        ivs = []
        for r in rows:
            s, e = sorted((int(r["gene_start"]), int(r["gene_end"])))
            ivs.append((s - WIN, e + WIN, r))
        ivs.sort(key=lambda x: x[0])
        cur = None
        groups = []
        for s, e, r in ivs:
            if cur and s <= cur[1]:
                cur = (cur[0], max(cur[1], e), cur[2] + [r])
            else:
                if cur:
                    groups.append(cur)
                cur = (s, e, [r])
        if cur:
            groups.append(cur)
        for s, e, members in groups:
            bi += 1
            bid = "blk%06d" % bi
            wrapped = False
            trunc = False
            if s < 1 or e > L:
                if circ and L:
                    wrapped = True
                else:
                    trunc = True
            fs, fe = max(1, s), min(L, e) if L else e
            path = os.path.join(seqdir, bid + ".fna")
            seq = None
            if not os.path.exists(path):
                seq, err = fetch_range(k, fs, fe)
                if seq is None:
                    failures.append((k, bid, err))
                else:
                    if wrapped and L:
                        if s < 1:
                            pre, _ = fetch_range(k, L + s + 1, L)
                            seq = (pre or "") + seq
                        if e > L:
                            post, _ = fetch_range(k, 1, e - L)
                            seq = seq + (post or "")
                    with open(path, "w", encoding="utf-8", newline="\n") as fh:
                        fh.write(">%s|%s|%d-%d\n" % (bid, k, fs, fe))
                        for j in range(0, len(seq), 70):
                            fh.write(seq[j:j + 70] + "\n")
            if seq is None and os.path.exists(path):
                seq = "".join(l.strip() for l in open(path, encoding="utf-8")
                              if not l.startswith(">"))
            h = hashlib.sha256((seq or "").encode()).hexdigest() if seq else ""
            blocks.append({"block_id": bid, "replicon_accession": k,
                           "replicon_length": L,
                           "topology": topo.get(k, "unknown"),
                           "block_start": fs, "block_end": fe,
                           "block_span_bp": len(seq or ""),
                           "wrapped_circular": "yes" if wrapped else "no",
                           "truncated": "yes" if trunc else "no",
                           "n_args_in_block": len(members),
                           "extracted_sequence_sha256": h})
            for r in members:
                s0, e0 = sorted((int(r["gene_start"]), int(r["gene_end"])))
                wrows.append({"block_id": bid,
                              "assembly_version": r["assembly_version"],
                              "replicon_accession": k, "replicon_length": L,
                              "topology": topo.get(k, "unknown"),
                              "determinant_name": r["determinant_name"],
                              "gene_family": r["gene_family"],
                              "gene_start": s0, "gene_end": e0,
                              "strand": r["strand"],
                              "window_start": s0 - WIN, "window_end": e0 + WIN,
                              "window_bp": WIN,
                              "wrapped_circular": "yes" if wrapped else "no",
                              "truncated": "yes" if trunc else "no",
                              "block_sequence_sha256": h})
        if bi % 2000 < len(groups):
            print("  blocks %d" % bi, flush=True)

    write(os.path.join(out, "shared_context_blocks.tsv"), blocks, list(blocks[0].keys()))
    write(os.path.join(out, "arg_neighbourhood_windows.tsv"), wrows, list(wrows[0].keys()))
    man = [{"window_set": "chromosomal_primary_10kb", "n_blocks": len(blocks),
            "n_arg_windows": len(wrows), "window_bp": WIN,
            "sensitivity_windows_bp": "5000;20000",
            "n_wrapped": sum(1 for b in blocks if b["wrapped_circular"] == "yes"),
            "n_truncated": sum(1 for b in blocks if b["truncated"] == "yes"),
            "total_extracted_bp": sum(b["block_span_bp"] for b in blocks),
            "n_failures": len(failures),
            "plasmid_neighbourhoods": "not retrieved separately; whole plasmid held"}]
    write(os.path.join(out, "arg_neighbourhood_manifest.tsv"), man, list(man[0].keys()))
    print("\nblocks %d | ARG windows %d | wrapped %d | truncated %d | failures %d"
          % (len(blocks), len(wrows), man[0]["n_wrapped"], man[0]["n_truncated"],
             len(failures)))
    print("extracted %.3f GB" % (man[0]["total_extracted_bp"] / 1e9))
    if failures:
        print("FAILURES:", failures[:5])


if __name__ == "__main__":
    main()
