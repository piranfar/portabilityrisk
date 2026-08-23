"""PR-CONTEXT Phase B -- BioSample epidemiological metadata for the primary cohort.

Runs on the LAPTOP because that is where the NCBI API key lives. The credential is read from
the environment, never printed, never written to any output file, and never sent to the
server. Only the resulting public metadata crosses the wire.

The extraction rule is deliberately conservative. A field is populated only when the
BioSample record explicitly documents it. Nothing is inferred: a host is not guessed from a
species name, a clinical status is not guessed from a hospital-sounding project title, and a
country is never refined into a city. Where the record is silent the field stays empty and
the genome counts against that field's coverage denominator.

BioSample-to-assembly linkage is preserved exactly, including the non-1:1 cases -- the cohort
has 6,288 genomes across 6,285 BioSamples, and those three shared BioSamples are kept and
flagged rather than collapsed.
"""
import argparse, collections, csv, hashlib, json, os, re, sys, time
import urllib.parse, urllib.request, datetime
import xml.etree.ElementTree as ET

VERSION = "pr_context_biosample_metadata_v1.0.0"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 200
_KEY = os.environ.get("NCBI_API_KEY")
_GAP = 0.11 if _KEY else 0.36

# Attribute harmonisation. Keys are lowercased attribute names as they appear in BioSample.
COUNTRY = {"geo_loc_name", "geographic location (country and/or sea)", "country",
           "geo loc name", "geographic location"}
COLL = {"collection_date", "collection date", "sample collection date"}
HOST = {"host", "host scientific name", "specific_host", "host_common_name"}
SOURCE = {"isolation_source", "isolation source", "isolate_source", "source_type",
          "sample_type", "env_medium", "environmental_medium", "host_tissue_sampled",
          "body_site", "host_body_site", "source"}

# Source-context classification. Every pattern below is a DOCUMENTED word in the record.
CLINICAL = re.compile(r"\b(clinical|patient|hospital|blood|urine|sputum|wound|csf|"
                      r"cerebrospinal|catheter|abscess|pus|bronch|tracheal|rectal swab|"
                      r"perianal|nosocomial|icu|bacteremia|septic|infection|isolate from "
                      r"patient|throat swab|body fluid|drain fluid|surgical)\b", re.I)
HUMAN = re.compile(r"\b(homo sapiens|human|patient)\b", re.I)
ANIMAL = re.compile(r"\b(bos taurus|sus scrofa|gallus|canis|felis|equus|ovis|animal|"
                    r"chicken|pig|swine|cattle|bovine|poultry|dog|cat|horse|sheep|"
                    r"turkey|duck|cow|calf|piglet|broiler|livestock|veterinary)\b", re.I)
FOOD = re.compile(r"\b(food|meat|retail|chicken meat|pork|beef|vegetable|salad|cheese|"
                  r"milk|seafood|shrimp|fish fillet|produce|lettuce|sprout)\b", re.I)
WASTE = re.compile(r"\b(wastewater|waste water|sewage|sewer|effluent|influent|"
                   r"treatment plant|wwtp|sludge)\b", re.I)
ENVIRON = re.compile(r"\b(soil|water|river|lake|sea|marine|sediment|environment|"
                     r"environmental|air|plant|rhizosphere|surface water|groundwater|"
                     r"hospital environment|sink|drain)\b", re.I)


def _post(ids, attempts=5):
    data = {"db": "biosample", "id": ",".join(ids), "retmode": "xml"}
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
                return r.read(), ""
        except Exception as e:
            last = repr(e)[:140]
            time.sleep(min(2 ** a, 30))
    return None, last


def date_precision(s):
    s = (s or "").strip()
    if not s:
        return "", ""
    if re.match(r"^\d{4}-\d{2}-\d{2}", s):
        return s[:10], "day"
    if re.match(r"^\d{4}-\d{2}$", s):
        return s, "month"
    if re.match(r"^\d{4}$", s):
        return s, "year"
    m = re.match(r"^(\d{4})", s)
    if m:
        return m.group(1), "year_from_freetext"
    return s, "unparsed"


def classify(source, host):
    """Documented-context classification. Silent record -> empty, never a guess."""
    blob = " ".join(x for x in (source, host) if x)
    if not blob.strip():
        return "", ""
    ctx = []
    if WASTE.search(blob):
        ctx.append("wastewater")
    if FOOD.search(blob):
        ctx.append("food")
    if ANIMAL.search(blob):
        ctx.append("animal")
    if HUMAN.search(blob) or CLINICAL.search(blob):
        ctx.append("human")
    if ENVIRON.search(blob) and "wastewater" not in ctx:
        ctx.append("environmental")
    clin = ""
    if CLINICAL.search(blob):
        clin = "clinical"
    elif ctx and "human" not in ctx:
        clin = "non_clinical"
    return ";".join(dict.fromkeys(ctx)), clin


COLS = ["biosample_accession", "assembly_versions", "n_assemblies_sharing_biosample",
        "bioproject_accession", "organism_harmonized", "genus", "country_raw", "country",
        "region_or_locality", "collection_date_raw", "collection_date", "date_precision",
        "host_raw", "host", "isolation_source_raw", "source_context",
        "clinical_status", "record_state"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")
    cache = os.path.join(a.root, "biosample_cache")
    os.makedirs(cache, exist_ok=True)
    print("%s | api key: %s (never printed, never leaves this machine)"
          % (VERSION, "present" if _KEY else "ABSENT"))

    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"]
    link = collections.defaultdict(set)
    meta = {}
    for r in occ:
        link[r["biosample_accession"]].add(r["assembly_version"])
        meta[r["biosample_accession"]] = (r["bioproject_accession"],
                                          r["organism_harmonized"], r["genus"])
    bss = sorted(x for x in link if x)
    print("unique BioSamples: %d | genomes: %d | shared-BioSample cases: %d"
          % (len(bss), len({r["assembly_version"] for r in occ}),
             sum(1 for b in bss if len(link[b]) > 1)))

    recs = {}
    for i in range(0, len(bss), BATCH):
        chunk = bss[i:i + BATCH]
        cp = os.path.join(cache, "b%05d.xml" % i)
        if os.path.exists(cp) and os.path.getsize(cp) > 0:
            raw = open(cp, "rb").read()
        else:
            raw, err = _post(chunk)
            if raw is None:
                print("  batch %d FAILED: %s" % (i, err))
                continue
            open(cp, "wb").write(raw)
        try:
            root = ET.fromstring(raw)
        except ET.ParseError as e:
            print("  batch %d parse error: %s" % (i, e))
            continue
        for bs in root.iter("BioSample"):
            acc = bs.get("accession") or ""
            attrs = {}
            for at in bs.iter("Attribute"):
                nm = (at.get("harmonized_name") or at.get("attribute_name") or "").lower()
                if nm and at.text:
                    attrs[nm] = at.text.strip()
            recs[acc] = attrs
        print("  %d/%d biosamples parsed" % (len(recs), len(bss)), flush=True)

    rows = []
    for b in bss:
        at = recs.get(b)
        bp, org, gen = meta[b]
        if at is None:
            rows.append({c: "" for c in COLS} | {
                "biosample_accession": b,
                "assembly_versions": ";".join(sorted(link[b])),
                "n_assemblies_sharing_biosample": len(link[b]),
                "bioproject_accession": bp, "organism_harmonized": org, "genus": gen,
                "record_state": "NOT_RETRIEVED"})
            continue
        craw = next((at[k] for k in at if k in COUNTRY), "")
        country, region = "", ""
        if craw:
            parts = [p.strip() for p in craw.split(":")]
            country = parts[0]
            region = parts[1] if len(parts) > 1 else ""
        draw = next((at[k] for k in at if k in COLL), "")
        d, prec = date_precision(draw)
        hraw = next((at[k] for k in at if k in HOST), "")
        sraw = next((at[k] for k in at if k in SOURCE), "")
        ctx, clin = classify(sraw, hraw)
        rows.append({
            "biosample_accession": b,
            "assembly_versions": ";".join(sorted(link[b])),
            "n_assemblies_sharing_biosample": len(link[b]),
            "bioproject_accession": bp, "organism_harmonized": org, "genus": gen,
            "country_raw": craw, "country": country, "region_or_locality": region,
            "collection_date_raw": draw, "collection_date": d, "date_precision": prec,
            "host_raw": hraw, "host": hraw, "isolation_source_raw": sraw,
            "source_context": ctx, "clinical_status": clin,
            "record_state": "RETRIEVED"})

    p = os.path.join(out, "biosample_metadata.tsv")
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in COLS) + "\n")

    n = len(rows)
    def cov(f):
        return sum(1 for r in rows if str(r.get(f, "")).strip())
    fields = [("record retrieved", "record_state"), ("country", "country"),
              ("region or locality", "region_or_locality"),
              ("collection date", "collection_date"), ("host", "host"),
              ("isolation source", "isolation_source_raw"),
              ("source context", "source_context"),
              ("clinical status", "clinical_status")]
    crows = [{"field": lab, "n_with_value": cov(f), "n_biosamples": n,
              "coverage_pct": round(100.0 * cov(f) / n, 2),
              "usable_for_stratified_analysis":
                  "yes" if cov(f) / n >= 0.50 else "no - below the 50% prespecified floor"}
             for lab, f in fields]
    with open(os.path.join(out, "biosample_field_coverage.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write("field\tn_with_value\tn_biosamples\tcoverage_pct\t"
                 "usable_for_stratified_analysis\n")
        for r in crows:
            fh.write("\t".join(str(r[k]) for k in
                               ("field", "n_with_value", "n_biosamples", "coverage_pct",
                                "usable_for_stratified_analysis")) + "\n")
    print("\n=== BioSample field coverage (denominator %d BioSamples) ===" % n)
    for r in crows:
        print("  %-22s %5d  %6.2f%%  %s" % (r["field"], r["n_with_value"],
                                            r["coverage_pct"],
                                            r["usable_for_stratified_analysis"]))
    print("\ntop countries:", dict(collections.Counter(
        r["country"] for r in rows if r["country"]).most_common(8)))
    print("source contexts:", dict(collections.Counter(
        r["source_context"] for r in rows if r["source_context"]).most_common(8)))
    print("clinical status:", dict(collections.Counter(
        r["clinical_status"] for r in rows if r["clinical_status"]).most_common()))
    print("biosample_metadata.tsv sha256:",
          hashlib.sha256(open(p, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
