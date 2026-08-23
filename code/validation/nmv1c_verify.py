"""NM-V1C independent blinding and integrity verification.

Imports nothing from the builder. Parses the delivered application as text, extracts its
embedded payload, and checks it against the frozen protocol, the manifest and the sealed key
without opening the key's label columns for any purpose other than set comparison.
"""
import argparse, collections, csv, hashlib, json, os, re, sys

VERSION = "nmv1c_verify_v1.1.0"
PROTO_SHA = "b2058877c0f7165f0ab7bf9be7adcb2cb3b5b4c8d2b81143d85ab03fe0f0c04c"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    D = a.dir
    SFX = a.suffix
    fails = []

    def C(name, cond, detail=""):
        print("  %-56s %s %s" % (name, "PASS" if cond else "*** FAIL ***", detail))
        if not cond:
            fails.append(name)

    app_p = os.path.join(D, "NMV1C_ADJUDICATION_APP%s.html" % SFX)
    man = json.load(open(os.path.join(D, "NMV1C_SAMPLE_MANIFEST%s.json" % SFX),
                     encoding="utf-8"))
    key = list(csv.DictReader(open(os.path.join(D, "NMV1C_UNBLINDING_KEY.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    txt = open(app_p, encoding="utf-8").read()
    print("%s\n  application %d bytes" % (VERSION, len(txt)))

    C("frozen protocol digest unchanged",
      sha(os.path.join(D, "NMV1C_FROZEN_PROTOCOL.json")) == PROTO_SHA)

    m = re.search(r"const CASES=(\[.*?\]);const REASONS=", txt, re.S)
    C("embedded payload parses", bool(m))
    cases = json.loads(m.group(1))
    toks = [c["token"] for c in cases]
    C("exactly 120 cases in the application", len(cases) == 120, len(cases))
    C("120 unique tokens", len(set(toks)) == 120)
    C("token sets identical: app vs manifest count", len(toks) == man["n"])
    C("token sets identical: app vs sealed key",
      set(toks) == {r["token"] for r in key}, "%d vs %d" % (len(set(toks)), len(key)))
    C("manifest all_token_hash matches key blocks",
      hashlib.sha256("|".join(sorted(r["block_id"] for r in key)).encode()).hexdigest()
      == man["all_token_hash"])

    POOL = json.load(open(a.pool, encoding="utf-8"))
    excl = set(POOL["excluded"])
    blocks = {r["block_id"] for r in key}
    C("zero overlap with every previously delivered package",
      not (blocks & excl), "%d prior blocks" % len(excl))
    C("all 120 blocks come from the unused population", blocks <= set(POOL["available"]))

    # ---- Method X / HMM leakage ----
    XTERMS = ["methodX", "Method X", "method_x", "hmm", "HMM", "transposase_marker",
              "mge_feature_inventory", "hmm_is", "hmm_integron", "phmm", "pHMM",
              "hmmsearch", "phmmer", "IS6_", "IS4_", "Phage_integrase", "new_343"]
    hits = [t for t in XTERMS if t in txt]
    C("zero Method X / HMM terms anywhere in the application", not hits, hits or "")
    payload_txt = json.dumps(cases)
    C("zero Method X terms in the embedded payload",
      not [t for t in XTERMS if t in payload_txt])
    methods = {f["m"] for c in cases for f in c["features"]}
    C("only Method Y and Z appear as feature sources", methods <= {"Y", "Z"}, sorted(methods))

    # ---- identity leakage ----
    ev = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(D, "nmv1_block_evidence_table.tsv"), encoding="utf-8"),
        delimiter="\t")}
    gt = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(D, "NMV1_RULE_BASED_GROUND_TRUTH.tsv"), encoding="utf-8"),
        delimiter="\t")}
    ident = set()
    for b in blocks:
        ident.update({b, ev[b]["species"], ev[b]["bioproject"], ev[b]["stratum"],
                      gt[b]["rule_based_label"]})
    ident.discard("")
    leak = [x for x in ident if x and x in txt]
    C("zero block id / species / BioProject / stratum / label leakage",
      not leak, "%d strings tested" % len(ident))
    C("no rule_id field exposed", "rule_id" not in txt)
    C("no machine_state field exposed", "machine_state" not in txt)
    C("no accession-like strings in payload",
      not re.search(r"(NZ_|NC_|GC[FA]_|PRJ)[A-Z0-9_.]+", payload_txt))

    # ---- instrument behaviour ----
    C("all decisions initially blank",
      not any(c.get("decision") or c.get("reason") for c in cases))
    for ch in ("MOBILE", "QUIESCENT", "NON_EVALUABLE"):
        C("choice present: %s" % ch, ch in txt)
    P = json.load(open(os.path.join(D, "NMV1C_FROZEN_PROTOCOL.json"), encoding="utf-8"))
    rc = P["instrument"]["reason_codes"]
    C("all %d reason codes present" % len(rc), all(r in txt for r in rc))
    C("reason enforced for NON_EVALUABLE",
      'd.decision!="NON_EVALUABLE"||d.reason' in txt.replace(" ", ""))
    C("keyboard shortcuts wired", 'e.key=="1"' in txt and 'ArrowRight' in txt)
    C("progress indicator present", "id='prog'" in txt or 'id="prog"' in txt)
    C("local autosave present", "localStorage" in txt)
    C("TSV export present", "NMV1C_ADJUDICATED_120.tsv" in txt)
    C("JSON export present", "NMV1C_ADJUDICATED_120.json" in txt)
    C("export carries only token/decision/reason/note",
      "token\\tdecision\\treason_code\\tnote" in txt)
    C("export cannot contain unblinded identifiers",
      not re.search(r"block_id|species|bioproject|rule_", txt.split("function exp(")[1]))
    C("offline: no external resource is requested",
      not re.search(r"(src|href)\s*=\s*[\"']https?://", txt))
    C("feature maps rendered for every case", txt.count("function svg(") == 1
      and all("features" in c for c in cases))

    # ---- stratum balance and weights ----
    st = collections.Counter(r["stratum"] for r in key)
    for s, v in P["selection"]["allocation"].items():
        C("stratum %s drawn as allocated" % s, st.get(s, 0) == v, "%d/%d" % (st.get(s, 0), v))
    C("inverse-probability weights retained in the key",
      all(float(r["weight"]) > 0 for r in key))
    ms = collections.Counter(r["machine_state"] for r in key)
    C("sampled states are MOBILE and QUIESCENT only",
      set(ms) == {"MOBILE", "QUIESCENT"}, dict(ms))

    # ---- start-up robustness (the R1 defect: a storage failure blanked the page) ----
    C("storage access is wrapped, never fatal",
      "function lget()" in txt and "function lset(" in txt)
    C("no unguarded storage call remains",
      len(re.findall(r"localStorage\.(?:get|set)Item", txt)) == 2
      and "function lget()" in txt)
    C("initial render is guarded against any exception",
      "try{render();}" in txt.replace("\n", "").replace(" ", ""))
    C("a visible banner reports lost autosave",
      "id='nostore'" in txt and "function banner()" in txt)
    C("progress can be restored from an export",
      "function imp(" in txt and "FileReader" in txt)
    _imp = txt.split("function imp(")
    C("restore accepts only the four exported fields",
      len(_imp) > 1 and not re.search(r"block_id|species|bioproject|stratum|rule_",
                                      _imp[1].split("try{render();}")[0]))

    print("\n  VERDICT: %s" % ("PASS - zero failures" if not fails
                               else "*** FAIL: %s ***" % fails))
    sys.exit(0 if not fails else 9)


if __name__ == "__main__":
    main()
