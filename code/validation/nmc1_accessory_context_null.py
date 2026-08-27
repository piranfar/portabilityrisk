"""Referee item 4.4: does the enrichment survive when the null keeps the neighbourhood?

The registered null relocates each chromosomal ARG occurrence uniformly across
its own chromosome. The paper concedes that this does not preserve accessory-
genome location: resistance genes sit in accessory regions, insertion sequences
accumulate in the same regions, and a null that can put an occurrence anywhere on
the chromosome will mostly put it in core sequence, where elements are sparse.
The enrichment could then be a statement about where accessory DNA is, not about
where resistance genes are within it.

Two arms, both registered here before any outcome is read.

ARM LOC - a LOCAL relocation null.
    Identical to the registered null in every respect except the relocation
    range: each occurrence moves uniformly within +/-R of its own position
    instead of anywhere on the chromosome. At R = 100 kb the null draws stay
    inside the occurrence's own neighbourhood - the same accessory island, the
    same element-dense region - so the accessory-location explanation is
    largely absorbed into the null itself. R = 50 kb is registered at the same
    time as a second radius, not chosen afterwards.

    One deliberate difference from C1-4: relocated intervals are NOT forced to
    be mutually non-overlapping. Inside a +/-50 kb window around a cluster of
    resistance genes the constraint is often unsatisfiable, and it does not
    enter this estimator anyway - each occurrence's distance depends only on its
    own relocated interval and the frozen element coordinates, never on where
    the other occurrences landed.

ARM ISL - a resistance-island stratification.
    Occurrences are single-linkage clustered along each replicon at a 10 kb gap.
    A cluster of >= 3 occurrences is treated as a resistance island; everything
    else is isolated. The enrichment is then recomputed within each stratum from
    the EXISTING null matrices, with no new permutation.

    This is a structural proxy, and the paper must say so: it is not a core
    versus accessory partition, and it does not identify AbaR or AbGRI islands
    by their backbone. It asks the narrower question the referee's second option
    asks - whether the effect is carried by multi-gene resistance islands - and
    it answers only that.

What is NOT done, stated rather than glossed: no core-genome alignment exists in
this project, so no true accessory-restricted null was run. The local null is the
closest available approximation and is labelled as an approximation.

ISEScan is never rerun. The element coordinate table is the frozen census output.
"""
import collections
import csv
import hashlib
import io
import json
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
csv.field_size_limit(1 << 30)

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", ".."))
NM = os.path.join(REPO, "docs", "nature_microbiology")
DATA = os.path.abspath(os.path.join(REPO, "..", "_nmbg_c1"))

AMEND = os.path.join(NM, "NM_C1_ACCESSORY_CONTEXT_AMENDMENT_009.json")
OUT = os.path.join(NM, "NM_C1_ACCESSORY_CONTEXT_RESULTS_V1.json")
MATS = os.path.join(DATA, "c1_4_matrices.npz")
ELEM = os.path.join(DATA, "genome_wide_is_elements.tsv")
MAN = os.path.join(NM, "NMBG_CHROMOSOME_ACCESSION_MANIFEST_FROZEN_V1.tsv")

HORIZON = 10000
LANDMARK = 1000
NPERM = 2000
SEED = 20260827
RADII = (50000, 100000)
ISLAND_GAP = 10000
ISLAND_MIN = 3

RATIO_FLOOR = 1.5          # reused from amendments 006-008, not chosen here
GROUPS = ("A. baumannii", "Klebsiella group", "P. aeruginosa", "Enterobacter group")

BODY = {
 "amendment": "NM-C1-009",
 "title": "Accessory-context arms for the background-normalised enrichment",
 "raised_by": "the uniform relocation null does not preserve accessory-genome location; "
              "resistance genes and insertion sequences may both concentrate in accessory "
              "regions, in which case a chromosome-wide null overstates the enrichment",
 "frozen_before_any_outcome_was_scored": True,
 "arms": {
   "LOC": {"name": "local relocation null",
           "rule": "each chromosomal occurrence is relocated uniformly within +/-R of its "
                   "own start instead of anywhere on its chromosome; interval length, "
                   "chromosome, genome, species and BioProject identity preserved; circular "
                   "chromosomes wrap; linear chromosomes exclude positions whose +/-10 kb "
                   "window would run past an end, exactly as in NM-C1-001",
           "radii_bp": list(RADII),
           "both_radii_registered_now": "R = 50 kb and R = 100 kb are fixed here together, "
                                        "so neither can be selected after seeing the other",
           "deliberate_difference_from_NM_C1_001":
               "relocated intervals are not forced to be mutually non-overlapping. Within a "
               "local window around a resistance-gene cluster the constraint is frequently "
               "unsatisfiable, and it does not enter this estimator: each occurrence's "
               "distance depends only on its own interval and the frozen element coordinates.",
           "permutations": NPERM, "seed": SEED,
           "per_chromosome_seed": "SEED XOR int(sha256(replicon_accession)[:16], 16)",
           "estimator": "unchanged - F(1 kb) under the NMIS complete-structural element "
                        "definition, right-censored at 10 kb",
           "isescan_rerun": False},
   "ISL": {"name": "resistance-island stratification",
           "rule": "occurrences on a replicon are single-linkage clustered at a %d bp gap; "
                   "a cluster of >= %d occurrences is an island, everything else is isolated"
                   % (ISLAND_GAP, ISLAND_MIN),
           "null": "the existing NM-C1-004 matrices, unchanged; no new permutation",
           "what_it_is_not": "not a core/accessory partition and not AbaR/AbGRI identification "
                             "by backbone. A structural proxy for multi-gene resistance "
                             "islands, reported as such."}},
 "scope_limit": "no core-genome alignment exists in this project, so no accessory-restricted "
                "null was run. ARM LOC is the closest available approximation and is reported "
                "as an approximation, not as an accessory-restricted null.",
 "decision_rule": {
   "thresholds_reused_from": "NM-C1-006, 007 and 008, deliberately, rather than chosen now",
   "accessory_location_explains_the_effect_if":
       "A. baumannii's enrichment fails to exclude its null 95%% interval under ARM LOC at "
       "R = 100 kb, or fails to exclude it among ISOLATED occurrences under ARM ISL, or the "
       "A. baumannii / Klebsiella contrast falls below %s in either" % RATIO_FLOOR,
   "reported_regardless": "every arm, every radius and every stratum, whatever they show",
   "no_threshold_may_move": True}
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def freeze():
    if os.path.exists(AMEND):
        return
    body = json.dumps(BODY, indent=1, ensure_ascii=False, sort_keys=True)
    d = dict(BODY)
    d["sha256_of_body"] = hashlib.sha256(body.encode("utf-8")).hexdigest()
    io.open(AMEND, "w", encoding="utf-8", newline="\n").write(
        json.dumps(d, indent=1, ensure_ascii=False, sort_keys=True) + "\n")
    print("frozen: %s  body %s" % (os.path.basename(AMEND), d["sha256_of_body"][:16]))


def verify():
    doc = json.load(io.open(AMEND, encoding="utf-8"))
    claimed = doc.pop("sha256_of_body")
    got = hashlib.sha256(json.dumps(doc, indent=1, ensure_ascii=False,
                                    sort_keys=True).encode("utf-8")).hexdigest()
    if got != claimed:
        raise SystemExit("REFUSING: amendment body changed since freeze")
    return doc, claimed


def is_complete(r):
    """NMIS complete_structural_IS, verbatim from NM-C1-004."""
    if r.get("type", "").strip() != "c":
        return False
    try:
        if not (float(r.get("orfLen") or 0) > 0):
            return False
        for k in ("orfBegin", "orfEnd"):
            if not str(r.get(k, "")).strip():
                return False
        if not (float(r.get("irLen") or 0) > 0):
            return False
        for k in ("start1", "end1", "start2", "end2"):
            if not str(r.get(k, "")).strip():
                return False
    except (TypeError, ValueError):
        return False
    return True


def chrom_seed(acc):
    h = int(hashlib.sha256(acc.encode("utf-8")).hexdigest()[:16], 16)
    return (SEED ^ h) & 0xFFFFFFFF


def group_of(sp):
    if sp == "Acinetobacter baumannii":
        return "A. baumannii"
    if sp.startswith("Klebsiella"):
        return "Klebsiella group"
    if sp == "Pseudomonas aeruginosa":
        return "P. aeruginosa"
    if sp.startswith("Enterobacter"):
        return "Enterobacter group"
    return None


def local_null_distances(ob, oe, eb, ee, span, circular, radius, rng, nperm):
    """(nperm, n) distances under relocation within +/-radius of each own start.

    Returns np.inf where no complete element lies fully inside the +/-10 kb
    window of the relocated interval - right-censored, never 10,000.
    """
    n = ob.size
    length = oe - ob
    lo = ob - radius
    hi = ob + radius
    if not circular:
        lo = np.maximum(lo, HORIZON)
        hi = np.minimum(hi, span - HORIZON - length)
        bad = hi <= lo
        if bad.any():                       # window larger than the usable span
            lo = np.where(bad, 0, lo)
            hi = np.where(bad, np.maximum(1, span - length), hi)
    out = np.empty((nperm, n), dtype=np.float64)
    if eb.size == 0:
        out[:] = np.inf
        return out
    # chunk over permutations so (chunk, n, m) stays small
    m = eb.size
    chunk = max(1, int(2e7 // max(1, n * m)))
    for a in range(0, nperm, chunk):
        b = min(nperm, a + chunk)
        s = rng.integers(lo, np.maximum(lo + 1, hi + 1), size=(b - a, n))
        if circular:
            s = np.mod(s, span)
        e = s + length
        ws = s[:, :, None] - HORIZON
        we = e[:, :, None] + HORIZON
        inside = (eb[None, None, :] >= ws) & (ee[None, None, :] <= we)
        gap = np.maximum(eb[None, None, :] - e[:, :, None],
                         s[:, :, None] - ee[None, None, :])
        gap = np.maximum(gap, 0)
        gap = np.where(inside, gap, np.inf)
        d = gap.min(axis=2)
        out[a:b] = np.where(d > HORIZON, np.inf, d)
    return out


def enrich(obs, null, mask):
    """F(1 kb) observed, expected and enrichment over the masked occurrences."""
    if mask.sum() == 0:
        return None
    o = float(((obs[mask] >= 0) & (obs[mask] <= LANDMARK)).mean())
    per = ((null[:, mask] >= 0) & (null[:, mask] <= LANDMARK)).mean(axis=1)
    e = float(per.mean())
    lo, hi = (float(x) for x in np.percentile(per, [2.5, 97.5]))
    return {"n": int(mask.sum()), "observed": round(o, 6), "expected": round(e, 6),
            "null_ci95": [round(lo, 6), round(hi, 6)],
            "excludes_null": bool(o < lo or o > hi),
            "enrichment": round(o / e, 4) if e > 0 else None}


def main():
    freeze()
    doc, claimed = verify()

    z = np.load(MATS, allow_pickle=False)
    ids, acc, species, assembly = z["ids"], z["acc"], z["species"], z["assembly"]
    obs = z["obs"].astype(np.float64)
    null0 = z["null"].astype(np.float64)
    # int16 matrices store right-censored occurrences as a sentinel; C1-4 wrote
    # -1 for "no qualifying element". Restore it as censored, not as distance 0.
    obs = np.where(obs < 0, np.inf, obs)
    null0 = np.where(null0 < 0, np.inf, null0)
    grp = np.array([group_of(s) for s in species])
    print("occurrences: %s   groups: %s"
          % ("{:,}".format(ids.size),
             ", ".join("%s %s" % (g, "{:,}".format(int((grp == g).sum()))) for g in GROUPS)))

    start = np.array([int(i.split("|")[-2]) for i in ids])
    end = np.array([int(i.split("|")[-1]) for i in ids])

    # ---- ARM ISL ---------------------------------------------------------
    island = np.zeros(ids.size, dtype=bool)
    by_rep = collections.defaultdict(list)
    for i in range(ids.size):
        by_rep[acc[i]].append(i)
    for _rep, idx in by_rep.items():
        idx = sorted(idx, key=lambda k: start[k])
        run = [idx[0]]
        for k in idx[1:]:
            if start[k] - end[run[-1]] <= ISLAND_GAP:
                run.append(k)
            else:
                if len(run) >= ISLAND_MIN:
                    island[run] = True
                run = [k]
        if len(run) >= ISLAND_MIN:
            island[run] = True
    print("island occurrences: %s of %s (%.1f%%)"
          % ("{:,}".format(int(island.sum())), "{:,}".format(ids.size),
             100 * island.mean()))

    isl = {}
    print("\nARM ISL - resistance-island stratification, existing null")
    print("  %-20s %-9s %8s %10s %10s %11s %s"
          % ("group", "stratum", "n", "observed", "expected", "enrichment", "excl null"))
    for g in GROUPS:
        isl[g] = {}
        for name, m in (("island", (grp == g) & island), ("isolated", (grp == g) & ~island)):
            r = enrich(obs, null0, m)
            isl[g][name] = r
            if r:
                print("  %-20s %-9s %8s %10.4f %10.4f %11.2f %s"
                      % (g, name, "{:,}".format(r["n"]), r["observed"], r["expected"],
                         r["enrichment"], "yes" if r["excludes_null"] else "NO"))
            else:
                print("  %-20s %-9s %8s   no occurrences" % (g, name, "0"))

    # ---- ARM LOC ---------------------------------------------------------
    man = {}
    with io.open(MAN, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            man[r["replicon_accession"]] = (int(r["replicon_length"]),
                                            r.get("topology", "").strip().lower() == "circular")
    el = collections.defaultdict(list)
    kept = 0
    with io.open(ELEM, encoding="utf-8", errors="replace", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if is_complete(r):
                el[r["replicon_accession"]].append((int(r["isBegin"]), int(r["isEnd"])))
                kept += 1
    print("\ncomplete structural elements loaded: %s on %s replicons"
          % ("{:,}".format(kept), "{:,}".format(len(el))))

    loc = {}
    for radius in RADII:
        nl = np.full((NPERM, ids.size), np.inf)
        missing_span = 0
        for rep, idx in by_rep.items():
            idx = np.array(sorted(idx))
            if rep not in man:
                missing_span += idx.size
                continue
            span, circular = man[rep]
            e = np.array(el.get(rep, []), dtype=np.int64)
            if e.size:
                eb, ee = e[:, 0], e[:, 1]
                if circular:      # replicate for wrap, as C1-4 does
                    eb = np.concatenate([eb - span, eb, eb + span])
                    ee = np.concatenate([ee - span, ee, ee + span])
            else:
                eb = ee = np.array([], dtype=np.int64)
            rng = np.random.default_rng(chrom_seed(rep))
            nl[:, idx] = local_null_distances(start[idx], end[idx], eb, ee,
                                              span, circular, radius, rng, NPERM)
        if missing_span:
            print("  WARNING: %d occurrences on replicons absent from the manifest"
                  % missing_span)
        res = {}
        print("\nARM LOC - local relocation null, R = %s kb" % "{:,}".format(radius // 1000))
        print("  %-20s %8s %10s %10s %11s %s"
              % ("group", "n", "observed", "expected", "enrichment", "excl null"))
        for g in GROUPS:
            r = enrich(obs, nl, grp == g)
            res[g] = r
            print("  %-20s %8s %10.4f %10.4f %11.2f %s"
                  % (g, "{:,}".format(r["n"]), r["observed"], r["expected"],
                     r["enrichment"], "yes" if r["excludes_null"] else "NO"))
        res["AB_over_Klebsiella"] = round(
            res["A. baumannii"]["enrichment"] / res["Klebsiella group"]["enrichment"], 3)
        print("  %-20s %48.2f" % ("A. baumannii / Klebsiella", res["AB_over_Klebsiella"]))
        loc["R%d" % radius] = res

    # ---- the registered rule --------------------------------------------
    ab_iso = isl["A. baumannii"]["isolated"]
    kl_iso = isl["Klebsiella group"]["isolated"]
    iso_ratio = round(ab_iso["enrichment"] / kl_iso["enrichment"], 3)
    r100 = loc["R100000"]
    fired = (not r100["A. baumannii"]["excludes_null"]
             or not ab_iso["excludes_null"]
             or r100["AB_over_Klebsiella"] < RATIO_FLOOR
             or iso_ratio < RATIO_FLOOR)
    print("\nA. baumannii / Klebsiella, isolated occurrences only: %.2f (floor %.1f)"
          % (iso_ratio, RATIO_FLOOR))
    print("frozen decision rule -> accessory location explains the effect: %s" % fired)

    rec = {"receipt": "NM-C1-009 accessory-context arms",
           "amendment_sha256_of_body": claimed, "amendment_sha256_of_file": sha(AMEND),
           "inputs": {"c1_4_matrices.npz": sha(MATS),
                      "genome_wide_is_elements.tsv": sha(ELEM),
                      "NMBG_CHROMOSOME_ACCESSION_MANIFEST_FROZEN_V1.tsv": sha(MAN)},
           "island_definition": {"gap_bp": ISLAND_GAP, "min_occurrences": ISLAND_MIN,
                                 "island_occurrences": int(island.sum()),
                                 "total_occurrences": int(ids.size)},
           "ARM_ISL": isl,
           "ARM_ISL_isolated_AB_over_Klebsiella": iso_ratio,
           "ARM_LOC": loc,
           "decision_rule": doc["decision_rule"],
           "accessory_location_explains_the_effect": bool(fired)}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("receipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
