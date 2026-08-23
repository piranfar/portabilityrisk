"""NMIS Phase 4 -- independent verification. Imports nothing from the scorer.

Re-derives census counts, structural predicates, coordinate conversion, containment, the crossing
category, denominators, estimates and contrasts with its own implementation, then compares to the
exported tables. Includes deterministic manual spot checks of all six required kinds.
"""
import argparse, collections, csv, hashlib, json, os, sys

VERSION = "nmis_verify_v1.0.0"
PROTO = "5438045a3b73d123347fcd60b2456779f050c2a8fae9cd016782b8b6168b03a3"
AMEND = "72ce50c2cfd512c1bc6cb74d29830f70474ec5f1a12056a9e744883246e67f80"
L = 10000
LAND = [1000, 2000, 5000, 10000]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def tsv(p):
    with open(p, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            if ln.strip():
                yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    D = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    NM = a.dir
    fails, n = [], 0

    def C(name, ok, det=""):
        nonlocal n
        n += 1
        print("  %-58s %s %s" % (name, "PASS" if ok else "*** FAIL ***", det))
        if not ok:
            fails.append(name)

    print("%s  (independent code path)\n" % VERSION)
    C("frozen protocol digest unchanged",
      sha(os.path.join(NM, "NMIS_FROZEN_PROTOCOL_V1.json")) == PROTO)
    C("amendment 001 digest unchanged",
      sha(os.path.join(NM, "NMIS_AMENDMENT_001_OCCURRENCE_WINDOW_CONTAINMENT.json")) == AMEND)
    R = json.load(open(os.path.join(NM, "NMIS_RESULT_RECEIPT_V1.json"), encoding="utf-8"))

    # ---- census integrity, recomputed ----
    st = list(tsv(os.path.join(a.census, "nmis_census_status.tsv")))
    ISD = os.path.join(a.census, "isescan")
    tsvs = sorted(f for f in os.listdir(ISD) if f.endswith(".tsv"))
    C("status rows = 21,955", len(st) == 21955, len(st))
    C("status block ids unique", len({r["block_id"] for r in st}) == 21955)
    C("tool failures = 0", sum(1 for r in st if r["status"] == "TOOL_FAILURE") == 0)
    C("result files = 5,563", len(tsvs) == 5563, len(tsvs))
    C("result files exactly match OK blocks",
      {f[:-4] for f in tsvs} == {r["block_id"] for r in st if r["status"] == "OK"})

    # ---- header-based parse, independent ----
    hdrsets, nel, types = set(), 0, collections.Counter()
    E = collections.defaultdict(list)
    for f in tsvs:
        with open(os.path.join(ISD, f), encoding="utf-8") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            hdrsets.add(tuple(head))
            for ln in fh:
                if not ln.strip():
                    continue
                v = ln.rstrip("\n").split("\t")
                rec = dict(zip(head, v))
                E[f[:-4]].append(rec)
                nel += 1
                types[rec["type"]] += 1
    C("single 24-column emitted header everywhere",
      len(hdrsets) == 1 and len(next(iter(hdrsets))) == 24)
    C("elements = 14,426", nel == 14426, nel)
    C("type c = 10,769, type p = 3,657",
      types["c"] == 10769 and types["p"] == 3657, dict(types))

    # ---- independent predicates and geometry ----
    blk = {r["block_id"]: r for r in tsv(os.path.join(D, "shared_context_blocks.tsv"))}
    win = list(tsv(os.path.join(D, "arg_neighbourhood_windows.tsv")))

    def I(v):
        v = (v or "").strip()
        if v in ("", "-", "NA", "nan"):
            return None
        try:
            return int(round(float(v)))
        except ValueError:
            return None

    def is_complete(e, blen):
        """Independent formulation of the frozen predicate."""
        if e["type"].strip() != "c":
            return False
        for k in ("orfBegin", "orfEnd", "start1", "end1", "start2", "end2"):
            if I(e[k]) is None:
                return False
        if not (I(e["orfLen"]) or 0) > 0:
            return False
        if not (I(e["irLen"]) or 0) > 0:
            return False
        if (I(e["isBegin"]) or 0) <= 1 or (I(e["isEnd"]) or 0) >= blen:
            return False          # boundary-touching stays partial
        return True

    def to_replicon(b, e):
        B = blk[b]
        off, RL = int(B["block_start"]), int(B["replicon_length"])
        s, t = I(e["isBegin"]), I(e["isEnd"])
        lo, hi = off + s - 1, off + t - 1
        if B["wrapped_circular"].strip().lower() in ("yes", "true", "1"):
            lo = (lo - 1) % RL + 1
            hi = (hi - 1) % RL + 1
            if hi < lo:
                return [(lo, RL), (1, hi)]
        return [(lo, hi)]

    def wsegs(w):
        s, t, RL = int(w["window_start"]), int(w["window_end"]), int(w["replicon_length"])
        return [(s, RL), (1, t)] if (w["wrapped_circular"].strip().lower() in
                                     ("yes", "true", "1") and t < s) else [(s, t)]

    def inside(ei, ws):
        for lo, hi in ei:
            if not any(lo >= s and hi <= t for s, t in ws):
                return False
        return True

    def sepr(a1, a2, b1, b2):
        return max(0, max(a1 - b2, b1 - a2))

    # ---- recompute endpoint states ----
    exported = {r["occurrence_id"]: r for r in
                tsv(os.path.join(NM, "nmis_occurrence_endpoints.tsv"))}
    C("exported endpoint table has 35,140 unique occurrences", len(exported) == 35140,
      len(exported))
    mism = 0
    dmism = 0
    cat = collections.Counter()
    ecat = collections.Counter()
    seen = set()
    spot = {"linear": None, "circular": None, "boundary": None, "crossing": None,
            "partial": None, "nooutput": None}
    for w in win:
        b = w["block_id"]
        oid = "|".join([w["assembly_version"], w["replicon_accession"], w["determinant_name"],
                        w["gene_start"], w["gene_end"]])
        B = blk[b]
        RL, blen = int(B["replicon_length"]), int(B["block_span_bp"])
        circ = w["topology"].strip().lower().startswith("circ")
        gs, ge = int(w["gene_start"]), int(w["gene_end"])
        ws = wsegs(w)
        best = None
        ncc = ncx = npt = 0
        for e in E.get(b, []):
            ei = to_replicon(b, e)
            comp = is_complete(e, blen)
            key = (b, e["isBegin"], e["isEnd"])
            if comp and inside(ei, ws):
                ncc += 1
                dd = min(min([sepr(gs, ge, lo, hi)] +
                             ([sepr(gs, ge, lo - RL, hi - RL), sepr(gs, ge, lo + RL, hi + RL)]
                              if circ else [])) for lo, hi in ei)
                if dd <= L:
                    best = dd if best is None else min(best, dd)
                if key not in seen:
                    ecat["complete_and_fully_contained"] += 1
                if spot["circular"] is None and B["wrapped_circular"].strip().lower() in \
                        ("yes", "true", "1"):
                    spot["circular"] = (b, e["isBegin"], e["isEnd"], ei, dd)
                if spot["linear"] is None and not circ:
                    spot["linear"] = (b, e["isBegin"], e["isEnd"], ei, dd)
            elif comp:
                ncx += 1
                if key not in seen:
                    ecat["complete_in_shared_block_but_crosses_occurrence_window_boundary"] += 1
                if spot["crossing"] is None:
                    spot["crossing"] = (b, e["isBegin"], e["isEnd"], ei, ws)
            else:
                npt += 1
                if key not in seen:
                    ecat["partial_or_boundary_limited"] += 1
                if spot["boundary"] is None and ((I(e["isBegin"]) or 0) <= 1 or
                                                 (I(e["isEnd"]) or 0) >= blen):
                    spot["boundary"] = (b, e["isBegin"], e["isEnd"], blen, e["type"])
                if spot["partial"] is None and e["type"].strip() == "p":
                    spot["partial"] = (b, e["isBegin"], e["isEnd"], e["type"])
            seen.add(key)
        state = ("complete_and_fully_contained" if best is not None else
                 "complete_in_shared_block_but_crosses_occurrence_window_boundary" if ncx else
                 "partial_or_boundary_limited" if npt else "no_structurally_resolved_IS")
        cat[state] += 1
        ex = exported[oid]
        if ex["endpoint_state"] != state:
            mism += 1
        exd = ex["nmis_dist_bp"]
        if (best is None) != (exd == ""):
            dmism += 1
        elif best is not None and int(exd) != best:
            dmism += 1
        if spot["nooutput"] is None and b not in E:
            spot["nooutput"] = (b, oid)
    C("endpoint states independently reproduced", mism == 0, "%d mismatches" % mism)
    C("structural distances independently reproduced", dmism == 0, "%d mismatches" % dmism)
    C("occurrence state counts match the receipt",
      cat == collections.Counter(R["occurrence_endpoint_states"]), dict(cat))
    C("element categories match the receipt",
      ecat == collections.Counter(R["element_categories"]), dict(ecat))
    C("element categories sum to 14,426", sum(ecat.values()) == 14426, sum(ecat.values()))

    # ---- denominators ----
    cls = collections.Counter(r["portability_class"] for r in exported.values())
    C("class A = 18,837 and class B = 16,303 unchanged",
      cls["A"] == 18837 and cls["B"] == 16303, dict(cls))
    C("blocks = 21,955", len({r["block_id"] for r in exported.values()}) == 21955)
    C("every occurrence has exactly one endpoint state",
      sum(cat.values()) == 35140)

    # ---- estimates recomputed ----
    pe = {r["group"]: r for r in tsv(os.path.join(NM, "nmis_primary_estimates.tsv"))}
    worst = 0.0
    for g in ("A. baumannii", "P. aeruginosa", "Klebsiella group"):
        sel = [r for r in exported.values() if r["frozen_group"] == g]
        W = sum(float(r["weight_block_balanced"]) for r in sel)
        for t in LAND:
            f = sum(float(r["weight_block_balanced"]) for r in sel
                    if r["nmis_censored"] == "0" and int(r["nmis_dist_bp"]) <= t) / W
            worst = max(worst, abs(f - float(pe[g]["F_%d" % t])))
        rmd = sum(float(r["weight_block_balanced"]) *
                  (int(r["nmis_dist_bp"]) if r["nmis_censored"] == "0" else L)
                  for r in sel) / W
        worst = max(worst, abs(rmd - float(pe[g]["restricted_mean_distance_bp"])) / L)
    C("all landmark and RMD estimates reproduced", worst < 1e-9, "max dev %.2e" % worst)

    con = list(tsv(os.path.join(NM, "nmis_species_contrasts.tsv")))
    C("contrast differences equal value_a minus value_b",
      all(abs((float(c["value_a"]) - float(c["value_b"])) - float(c["difference"])) < 1e-12
          for c in con))
    C("Holm value present on every contrast row", all(c["p_holm"] != "" for c in con))
    C("gate verdict follows from the 1 kb and 2 kb cells",
      (R["gate_verdict"] == "SUCCESS") ==
      all(float(c["difference"]) > 0 and float(c["ci_lo"]) > 0 and float(c["p_holm"]) <= 0.05
          for c in con if int(c["landmark_bp"]) in (1000, 2000)))

    # ---- row uniqueness and digests ----
    for f in ("nmis_occurrence_endpoints.tsv", "nmis_primary_estimates.tsv",
              "nmis_species_contrasts.tsv", "nmis_sensitivity_results.tsv",
              "nmis_vs_nmdist_comparison.tsv"):
        rows = list(tsv(os.path.join(NM, f)))
        key = [tuple(sorted(r.items())) for r in rows]
        C("row uniqueness: %s" % f[:44], len(set(key)) == len(key), "%d rows" % len(rows))
    C("receipt output digests match files on disk",
      all(sha(os.path.join(NM, k)) == v for k, v in R["output_digests"].items()))

    # ---- deterministic manual spot checks ----
    print("\n  DETERMINISTIC SPOT CHECKS")
    C("spot: ordinary linear containment", spot["linear"] is not None,
      str(spot["linear"])[:78] if spot["linear"] else "")
    C("spot: circular wrapped containment", spot["circular"] is not None,
      str(spot["circular"])[:78] if spot["circular"] else "")
    C("spot: boundary-touching element stays partial", spot["boundary"] is not None,
      str(spot["boundary"])[:78] if spot["boundary"] else "")
    C("spot: complete element crossing the occurrence window", spot["crossing"] is not None,
      str(spot["crossing"])[:70] if spot["crossing"] else "")
    C("spot: partial element", spot["partial"] is not None,
      str(spot["partial"])[:78] if spot["partial"] else "")
    C("spot: block with no ISEScan output", spot["nooutput"] is not None,
      str(spot["nooutput"])[:78] if spot["nooutput"] else "")

    print("\n  checks run: %d   disagreements: %d" % (n, len(fails)))
    print("  VERDICT: %s" % ("PASS - zero disagreements" if not fails
                             else "*** FAIL: %s ***" % fails))
    sys.exit(0 if not fails else 9)


if __name__ == "__main__":
    main()
