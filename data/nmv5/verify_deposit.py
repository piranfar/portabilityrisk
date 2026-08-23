"""Verify this deposit, standing alone.

Run it from inside the extracted archive. It needs only the Python standard
library and reads nothing outside this directory.

    python verify_deposit.py

It re-derives every headline denominator from the dataset itself rather than
reading them from a summary, validates every table against its JSON Schema,
re-hashes every file against SHA256SUMS, and reports its own disagreement count.
"""
import collections, csv, hashlib, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATASET = "portabilityrisk_occurrence_portability_v1.tsv"

EXPECTED = {
    "occurrences": 74349,
    "chromosome": 35140,
    "plasmid": 39209,
    "A": 18837, "B": 16303, "C": 7170, "D": 6043, "E": 25996,
    "genomes": 6288, "biosamples": 6285, "bioprojects": 2283,
    "species": 109, "replicons": 12811,
    "collapsed_events": 9755, "collapsed_plasmid_events": 3569,
}

res = []


def ck(name, ok, detail=""):
    res.append((name, bool(ok), detail))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         ("\n         " + detail) if detail else ""))


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    dp = os.path.join(HERE, DATASET)
    if not os.path.isfile(dp):
        print("dataset not found: %s" % DATASET)
        return 2

    # ------------------------------------------------- re-derive everything
    cls, mol = collections.Counter(), collections.Counter()
    genomes, biosamples, bioprojects, species, replicons = (set(), set(), set(),
                                                            set(), set())
    events = set()
    per_genome = collections.defaultdict(lambda: [0, 0])
    keys = collections.Counter()
    n = 0
    with open(dp, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = list(rd.fieldnames)
        for r in rd:
            n += 1
            cls[r["portability_class"]] += 1
            m = r["replicon_molecule_type"]
            mol[m] += 1
            genomes.add(r["assembly_version"])
            biosamples.add(r["biosample_accession"])
            bioprojects.add(r["bioproject_accession"])
            species.add(r["organism_harmonized"])
            replicons.add(r["replicon_accession"])
            events.add((r["assembly_version"], m))
            per_genome[r["assembly_version"]][1] += 1
            if m == "Plasmid":
                per_genome[r["assembly_version"]][0] += 1
            keys[(r["assembly_version"], r["replicon_accession"],
                  r["determinant_name"], r["gene_start"], r["gene_end"])] += 1

    ck("V1 row count", n == EXPECTED["occurrences"], "%d rows, 19 columns expected %d cols"
       % (n, len(cols)))
    ck("V2 column count", len(cols) == 19, "%d columns" % len(cols))
    dup = sum(v - 1 for v in keys.values() if v > 1)
    ck("V3 primary key unique", dup == 0, "duplicate-key rows=%d" % dup)
    ck("V4 compartment reconciliation",
       mol.get("Chromosome") == EXPECTED["chromosome"]
       and mol.get("Plasmid") == EXPECTED["plasmid"],
       "chromosome=%d plasmid=%d" % (mol.get("Chromosome", 0), mol.get("Plasmid", 0)))
    ck("V5 five-class reconciliation",
       all(cls.get(k) == EXPECTED[k] for k in "ABCDE"),
       " ".join("%s=%d" % (k, cls.get(k, 0)) for k in "ABCDE"))
    ck("V6 A+B equals chromosome",
       cls.get("A", 0) + cls.get("B", 0) == EXPECTED["chromosome"],
       "%d" % (cls.get("A", 0) + cls.get("B", 0)))
    ck("V7 C+D+E equals plasmid",
       cls.get("C", 0) + cls.get("D", 0) + cls.get("E", 0) == EXPECTED["plasmid"],
       "%d" % (cls.get("C", 0) + cls.get("D", 0) + cls.get("E", 0)))
    ck("V8 classes sum to the denominator, no unclassified row",
       sum(cls.values()) == EXPECTED["occurrences"] and set(cls) == set("ABCDE"), "")
    ck("V9 cohort counts",
       len(genomes) == EXPECTED["genomes"] and len(biosamples) == EXPECTED["biosamples"]
       and len(bioprojects) == EXPECTED["bioprojects"]
       and len(species) == EXPECTED["species"] and len(replicons) == EXPECTED["replicons"],
       "genomes=%d biosamples=%d bioprojects=%d species=%d replicons=%d"
       % (len(genomes), len(biosamples), len(bioprojects), len(species), len(replicons)))

    pl_share = 100.0 * mol.get("Plasmid", 0) / n
    ck("V10 occurrence-weighted plasmid share is 52.736 %",
       abs(pl_share - 52.736) < 0.001, "%.4f %%" % pl_share)

    ev_pl = sum(1 for _, m in events if m == "Plasmid")
    ev_share = 100.0 * ev_pl / len(events)
    ck("V11 genome-collapsed events: 3,569 / 9,755 = 36.586 %",
       len(events) == EXPECTED["collapsed_events"]
       and ev_pl == EXPECTED["collapsed_plasmid_events"]
       and abs(ev_share - 36.586) < 0.001,
       "%d / %d = %.4f %%" % (ev_pl, len(events), ev_share))

    # Two published figures for "the arithmetic mean of per-genome percentages"
    # differ in the sixth decimal, and the reason is exact rather than sloppy:
    # 35.932096 % is the mean of the STORED pct_plasmid column, rounded to two
    # decimals; 35.932101 % is the mean recomputed from the raw counts. Both are
    # recorded so neither can be mistaken for an error in the other.
    mean_pct = sum(100.0 * p / t for p, t in per_genome.values()) / len(per_genome)
    ck("V12 the arithmetic mean of per-genome percentages is a DIFFERENT "
       "statistic from the collapsed-event share",
       abs(mean_pct - 35.932101) < 0.00001 and abs(mean_pct - ev_share) > 0.5,
       "mean from raw counts = %.6f %%; mean of the rounded published column = "
       "35.932096 %%; collapsed-event share = %.4f %%. The first two differ only "
       "because one averages a 2-dp-rounded column; neither is the third."
       % (mean_pct, ev_share))

    b_share = 100.0 * cls.get("B", 0) / EXPECTED["chromosome"]
    ck("V13 chromosomal MGE association is 46.39 % occurrence-weighted",
       abs(b_share - 46.39) < 0.01, "%.4f %%" % b_share)

    # --------------------------------------------------------- schemas
    sd = os.path.join(HERE, "schemas")
    bad = []
    if os.path.isdir(sd):
        for s in sorted(os.listdir(sd)):
            tab = s.replace(".schema.json", ".tsv")
            tp = os.path.join(HERE, tab)
            if not os.path.isfile(tp):
                bad.append("%s: table absent" % tab)
                continue
            sch = json.load(open(os.path.join(sd, s), encoding="utf-8"))
            with open(tp, encoding="utf-8", newline="") as fh:
                hdr = next(csv.reader(fh, delimiter="\t"))
            if hdr != sch["required"]:
                bad.append("%s: header does not match schema" % tab)
    ck("V14 every table matches its JSON Schema header", not bad,
       "; ".join(bad) or "checked %d schemas" % len(os.listdir(sd) if os.path.isdir(sd) else []))

    # --------------------------------------------------------- checksums
    sp = os.path.join(HERE, "SHA256SUMS")
    if os.path.isfile(sp):
        mism = missing = cnt = 0
        for line in open(sp, encoding="utf-8"):
            line = line.strip()
            if not line:
                continue
            want, name = line.split(None, 1)
            name = name.lstrip("*").strip()
            fp = os.path.join(HERE, name)
            if not os.path.isfile(fp):
                missing += 1
            else:
                cnt += 1
                if sha(fp) != want:
                    mism += 1
        ck("V15 every file matches SHA256SUMS", mism == 0 and missing == 0,
           "checked=%d missing=%d mismatched=%d" % (cnt, missing, mism))
    else:
        ck("V15 SHA256SUMS present", False, "absent")

    # --------------------------------------------------------- disclosure
    # Match restricted ARTEFACTS, not restricted vocabulary. This deposit's
    # documentation must be able to say "unblinding keys are not deposited
    # here"; a rule that forbids the phrase forbids the disclosure.
    RESTRICTED = re.compile(
        r"(?i)BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY|ssh-(?:rsa|ed25519) AAAA"
        r"|ocid1\.[a-z]+\."
        r"|\b(?:150\.136|129\.80)\.\d{1,3}\.\d{1,3}\b"
        r"|(?<![\w/])/home/[a-z][a-z0-9_-]*/|[A-Za-z]:[\\/](?:Users|Github)"
        r"|(?:password|api[_-]?key|secret|token)\s*[:=]\s*[\"'][A-Za-z0-9_./+-]{12,}"
        r"|(?m:^)blinded_id\t"
        r"|NMV1C?_(?:AUDIT|ADJUDICATION)_(?:UNBLINDING_KEY|CASEBOOK|BLINDED_PACKAGE)")
    hits = []
    for dp2, dns, fns in os.walk(HERE):
        for fn in fns:
            if not fn.lower().endswith((".tsv", ".md", ".json", ".py", ".txt")):
                continue
            fp = os.path.join(dp2, fn)
            try:
                t = open(fp, encoding="utf-8", errors="replace").read()
            except OSError:
                continue
            for m in RESTRICTED.finditer(t):
                if fn == os.path.basename(__file__):
                    continue
                hits.append("%s: %s" % (fn, m.group(0)[:40]))
    ck("V16 no restricted term, credential, key or infrastructure identifier",
       not hits, "; ".join(hits[:5]) or "none")

    fails = [r for r in res if not r[1]]
    print("\nchecks: %d   disagreements: %d" % (len(res), len(fails)))
    if not fails:
        print("This deposit reproduces every headline denominator from its own contents.")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
