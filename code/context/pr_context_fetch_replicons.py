"""PR-CONTEXT step 2 -- retrieve the documented replicon inventory for all 7,216 genomes.

METADATA ONLY. Not one nucleotide is downloaded here. NCBI's sequence report states, per
assembly, every assembled molecule: its RefSeq accession, its GenBank accession, its length,
and -- the field the whole analysis turns on -- assigned_molecule_location_type, which is
NCBI's own designation of Chromosome or Plasmid.

That designation is documented evidence about a closed replicon in a complete genome. It is
not a prediction, and nothing here predicts anything.

Rate-limited to NCBI's published limits. Every payload is cached to disk so a rerun costs
nothing and so the exact bytes behind every row survive.
"""
import argparse, collections, csv, hashlib, json, os, queue, sys, threading, time
import datetime, urllib.error, urllib.request

VERSION = "pr_context_fetch_replicons_v1.0.0"
EP = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/%s/sequence_reports"

_KEY = os.environ.get("NCBI_API_KEY")
_SLOTS = threading.Semaphore(8 if _KEY else 3)
_GAP = 0.12 if _KEY else 0.36
_last = [0.0]
_pace = threading.Lock()
_lock = threading.Lock()


def _paced():
    with _pace:
        d = _GAP - (time.time() - _last[0])
        if d > 0:
            time.sleep(d)
        _last[0] = time.time()


def fetch(acc, cache, attempts=5):
    cp = os.path.join(cache, acc + ".json")
    if os.path.exists(cp) and os.path.getsize(cp) > 0:
        return json.load(open(cp, encoding="utf-8")), "cached", ""
    url = EP % acc
    if _KEY:
        url += "?api_key=" + _KEY
    last = ""
    for a in range(attempts):
        try:
            with _SLOTS:
                _paced()
                req = urllib.request.Request(
                    url, headers={"User-Agent": "portabilityrisk-context/1.0"})
                with urllib.request.urlopen(req, timeout=180) as r:
                    raw = r.read()
            tmp = cp + ".tmp"
            open(tmp, "wb").write(raw)
            os.replace(tmp, cp)
            return json.loads(raw), "fetched", ""
        except urllib.error.HTTPError as e:
            last = "HTTP %s" % e.code
            if e.code not in (429, 500, 502, 503, 504) or a == attempts - 1:
                return None, "failed", last
            time.sleep(min(2 ** a, 30))
        except Exception as e:
            last = repr(e)[:120]
            if a == attempts - 1:
                return None, "failed", last
            time.sleep(min(2 ** a, 30))
    return None, "failed", last


COLS = ["assembly_accession", "sequence_accession", "refseq_accession",
        "genbank_accession", "replicon_name", "replicon_length",
        "replicon_molecule_type", "assembly_unit", "role", "gc_percent",
        "sequence_name", "fetch_state"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--amrfinder-dir", required=True)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")
    cache = os.path.join(a.root, "seqreport_cache")
    for d in (out, cache):
        os.makedirs(d, exist_ok=True)

    accs = sorted(f[:-len(".tsv.run.json")]
                  for f in os.listdir(a.amrfinder_dir) if f.endswith(".tsv.run.json"))
    print("%s | assemblies to resolve: %d | api_key: %s"
          % (VERSION, len(accs), "yes" if _KEY else "no"))
    if len(accs) != 7216:
        print("REFUSING: expected 7,216 processed assemblies, found %d" % len(accs))
        sys.exit(1)

    q = queue.Queue()
    for x in accs:
        q.put(x)
    rows, failures, done, states = [], [], [0], collections.Counter()

    def worker():
        while True:
            try:
                acc = q.get_nowait()
            except queue.Empty:
                return
            d, st, err = fetch(acc, cache)
            recs = (d.get("reports") if isinstance(d, dict) else None) or []
            local = []
            if d is None:
                local.append({"assembly_accession": acc, "fetch_state": "FETCH_FAILED:" + err})
            elif not recs:
                local.append({"assembly_accession": acc, "fetch_state": "NO_SEQUENCE_REPORT"})
            else:
                for r in recs:
                    rs = r.get("refseq_accession") or ""
                    gb = r.get("genbank_accession") or ""
                    local.append({
                        "assembly_accession": acc,
                        # the identifier AMRFinderPlus actually printed is the one in the
                        # FASTA we fed it: RefSeq for GCF_, GenBank for GCA_
                        "sequence_accession": (rs if acc.startswith("GCF_") else gb) or rs or gb,
                        "refseq_accession": rs, "genbank_accession": gb,
                        "replicon_name": r.get("chr_name") or "",
                        "replicon_length": r.get("length") or "",
                        "replicon_molecule_type":
                            r.get("assigned_molecule_location_type") or "",
                        "assembly_unit": r.get("assembly_unit") or "",
                        "role": r.get("role") or "",
                        "gc_percent": r.get("gc_percent") or "",
                        "sequence_name": r.get("sequence_name") or "",
                        "fetch_state": st})
            with _lock:
                rows.extend(local)
                states[st] += 1
                if d is None:
                    failures.append((acc, err))
                done[0] += 1
                if done[0] % 500 == 0 or done[0] == len(accs):
                    print("  %5d/%d  replicon rows %d  failures %d"
                          % (done[0], len(accs), len(rows), len(failures)), flush=True)

    ts = [threading.Thread(target=worker, daemon=True) for _ in range(a.workers)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()

    rows.sort(key=lambda r: (r["assembly_accession"], r.get("sequence_accession", "")))
    p = os.path.join(out, "replicon_inventory.tsv")
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in COLS) + "\n")

    mt = collections.Counter(r.get("replicon_molecule_type", "") for r in rows)
    print("\nreplicon rows: %d over %d assemblies" % (len(rows), len(accs)))
    print("fetch states:", dict(states))
    print("molecule types:")
    for k, v in mt.most_common():
        print("  %-24s %6d" % (k or "(blank)", v))
    print("assemblies with zero replicon records: %d"
          % len({r["assembly_accession"] for r in rows if not r.get("sequence_accession")}))
    if failures:
        print("FAILURES (%d): %s" % (len(failures), failures[:10]))
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    json.dump({"generated_utc": datetime.datetime.now(datetime.timezone.utc)
               .strftime("%Y-%m-%dT%H:%M:%SZ"), "builder": VERSION,
               "assemblies": len(accs), "replicon_rows": len(rows),
               "fetch_states": dict(states), "molecule_types": dict(mt),
               "failures": len(failures), "replicon_inventory_sha256": h,
               "statement": "Metadata only. No nucleotide sequence was downloaded."},
              open(os.path.join(out, "replicon_inventory_receipt.json"), "w",
                   encoding="utf-8", newline="\n"), indent=2)
    print("replicon_inventory.tsv sha256 %s" % h)


if __name__ == "__main__":
    main()
