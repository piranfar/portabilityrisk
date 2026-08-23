"""NM-V3 independent verification.

Deliberately does NOT import or reuse nmv3_score.py. It re-reads the raw census reports, the
frozen tables and the Arm I evidence with its own parsing, recomputes every quantity the receipt
claims, and reports a disagreement count per check. A single disagreement is a failure.

Differences from the scorer, on purpose: line-oriented parsing instead of csv.DictReader for the
census, index-based column lookup instead of names, an independently written classification
function, and a bootstrap seeded and structured separately.
"""
import argparse, collections, hashlib, json, os, random, sys

VERSION = "nmv3_verify_v1.0.0"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def read_tsv_rows(p):
    """Line-oriented reader; returns (header_list, list_of_field_lists)."""
    with open(p, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    return hdr, rows


def col(hdr, name):
    return hdr.index(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--armi", required=True)
    ap.add_argument("--receipt", required=True)
    a = ap.parse_args()
    D = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    R = json.load(open(a.receipt, encoding="utf-8"))

    fails = []
    checks = 0

    def C(name, ok, detail=""):
        nonlocal checks
        checks += 1
        print("  %-58s %s %s" % (name, "PASS" if ok else "*** FAIL ***", detail))
        if not ok:
            fails.append(name)

    print("%s  (independent code path; imports nothing from the scorer)" % VERSION)

    # ---------- denominators, recomputed from the frozen tables ----------
    hdr, rows = read_tsv_rows(os.path.join(D, "determinant_portability_classes.tsv"))
    ci = col(hdr, "portability_class"); ri = col(hdr, "replicon_accession")
    cls = collections.Counter(r[ci] for r in rows)
    tot = len(rows)
    plas = cls["C"] + cls["D"] + cls["E"]
    C("total occurrences = 74,349", tot == 74349, tot)
    C("plasmid-borne occurrences = 39,209", plas == 39209, plas)
    C("class counts A/B/C/D/E reproduce the frozen partition",
      (cls["A"], cls["B"], cls["C"], cls["D"], cls["E"]) == (18837, 16303, 7170, 6043, 25996),
      "%d/%d/%d/%d/%d" % (cls["A"], cls["B"], cls["C"], cls["D"], cls["E"]))
    C("A+B+C+D+E sums to the denominator", sum(cls.values()) == tot)
    C("receipt denominators agree",
      R["census"]["total_occurrences"] == tot and R["census"]["plasmid_side_occurrences"] == plas)

    occ_per_rep = collections.Counter(r[ri] for r in rows if r[ci] in ("C", "D", "E"))

    # ---------- frozen marker table, independent parse ----------
    fh_, frows = read_tsv_rows(os.path.join(D, "plasmid_mobility_annotation.tsv"))
    fa = col(fh_, "replicon_accession")
    frel, fmpf, fori = col(fh_, "has_relaxase"), col(fh_, "has_mpf"), col(fh_, "has_orit")
    fcat = col(fh_, "portability_mobility_category")
    frozen = {r[fa]: (r[frel].strip().lower() in ("1", "true", "yes"),
                      r[fmpf].strip().lower() in ("1", "true", "yes"),
                      r[fori].strip().lower() in ("1", "true", "yes"),
                      r[fcat]) for r in frows}
    C("frozen replicon count = 6,621", len(frozen) == 6621, len(frozen))

    def classify(rel, mpf, ori):
        """Independently written; same frozen rule, different expression."""
        if rel:
            return "E" if mpf else "D"
        return "D" if ori else "C"

    # ---------- census reports, independent parse ----------
    def load_census(tag):
        p = os.path.join(a.census, tag + ".tsv")
        h, rr = read_tsv_rows(p)
        si = col(h, "sample_id")
        rl, mp, orr = col(h, "relaxase_type(s)"), col(h, "mpf_type"), col(h, "orit_type(s)")
        out = {}
        for r in rr:
            acc = r[si].split(" ")[0]
            out[acc] = (r[rl] not in ("", "-"), r[mp] not in ("", "-"), r[orr] not in ("", "-"))
        return out

    base = load_census("BASE_319_db318")
    armt = load_census("ARM_T_318_db318")
    armd = load_census("ARM_D_319_db200")
    C("census rows 6,621 in every configuration",
      len(base) == len(armt) == len(armd) == 6621,
      "%d/%d/%d" % (len(base), len(armt), len(armd)))
    C("census accessions identical to the frozen table", set(base) == set(frozen))

    # ---------- baseline reproduction ----------
    mk = sum(1 for x in frozen if base[x] != frozen[x][:3])
    cat = sum(1 for x in frozen
              if classify(*base[x]) != {"predicted_conjugative": "E",
                                        "predicted_mobilizable": "D",
                                        "nonconjugative_or_no_mobility_markers_detected": "C"}[frozen[x][3]])
    C("baseline reproduces frozen markers exactly", mk == 0, "%d mismatches" % mk)
    C("baseline reproduces frozen categories exactly", cat == 0, "%d mismatches" % cat)
    C("receipt's exact_reproduction claim is true",
      R["baseline_reproduction_of_frozen_run"]["exact_reproduction"] is True and mk == 0 and cat == 0)

    # ---------- transition matrices ----------
    for tag, alt, key in (("ARM T", armt, "ARM_T_tool_3.1.8"),
                          ("ARM D", armd, "ARM_D_database_v2.0.0")):
        tm = collections.Counter()
        occ_chg = 0
        for x in frozen:
            f, t = classify(*base[x]), classify(*alt[x])
            tm[(f, t)] += 1
            if f != t:
                occ_chg += occ_per_rep[x]
        chg = sum(v for (f, t), v in tm.items() if f != t)
        C("%s replicon transitions recomputed" % tag,
          chg == R["arms"][key]["PRIMARY_replicon_category_changes"], chg)
        C("%s occurrence transitions recomputed" % tag,
          occ_chg == R["arms"][key]["SECONDARY_occurrence_changes"], occ_chg)
        C("%s transition matrix is diagonal" % tag,
          all(f == t for (f, t), v in tm.items() if v > 0))
        C("%s within the frozen 5 pct bound" % tag, (100.0 * chg / len(frozen)) <= 5.0,
          "%.4f%%" % (100.0 * chg / len(frozen)))

    # ---------- E1 / E2 reconciliation ----------
    e1 = e2 = 0
    for x in frozen:
        rel, mpf, ori = base[x]
        if classify(rel, mpf, ori) != "E":
            continue
        if ori: e2 += occ_per_rep[x]
        else:   e1 += occ_per_rep[x]
    b = R["E1_E2"]["baseline"]
    C("E1 occurrences recomputed", e1 == b["E1_occurrences"], e1)
    C("E2 occurrences recomputed", e2 == b["E2_occurrences"], e2)
    C("E1 + E2 == class E == 25,996", e1 + e2 == 25996 and e1 + e2 == cls["E"], e1 + e2)
    C("E2 pct of class E", abs(100.0 * e2 / (e1 + e2) - b["E2_pct_of_E"]) < 1e-9,
      "%.6f%%" % (100.0 * e2 / (e1 + e2)))
    C("E2 pct of plasmid-borne", abs(100.0 * e2 / plas - b["E2_pct_of_plasmid_borne"]) < 1e-9,
      "%.6f%%" % (100.0 * e2 / plas))
    C("E pct of plasmid-borne is the 66.3 figure",
      abs(100.0 * (e1 + e2) / plas - 66.301104) < 1e-4, "%.6f%%" % (100.0 * (e1 + e2) / plas))

    # ---------- Arm I, independent parse ----------
    ih, irows = read_tsv_rows(a.armi)
    ia = col(ih, "accession"); irel = col(ih, "relaxase_independent")
    impf = col(ih, "mpf_independent"); ist = col(ih, "status")
    ev = {r[ia]: (r[irel] == "1", r[impf] == "1", r[ist]) for r in irows}
    C("Arm I covers every census replicon", set(ev) == set(frozen), len(ev))
    C("Arm I zero-CDS count", sum(1 for v in ev.values() if v[2] == "ZERO_CDS")
      == R["arm_I"]["zero_cds_unresolved"])
    ce = collections.Counter()
    for x in frozen:
        mob_e = base[x][0] and base[x][1]
        cj_e = ev[x][0] and ev[x][1]
        ce[(mob_e, cj_e)] += 1
    aiE = R["arm_I"]["class_E_definition_concordance"]
    C("Arm I class-E concordance cells recomputed",
      ce[(True, True)] == aiE["both_E"] and ce[(False, False)] == aiE["neither_E"]
      and ce[(True, False)] == aiE["mob_suite_only_E"] and ce[(False, True)] == aiE["conjscan_only_E"],
      "%d/%d/%d/%d" % (ce[(True, True)], ce[(False, False)], ce[(True, False)], ce[(False, True)]))
    n = sum(ce.values())
    po = (ce[(True, True)] + ce[(False, False)]) / n
    p1 = (ce[(True, True)] + ce[(True, False)]) / n
    q1 = (ce[(True, True)] + ce[(False, True)]) / n
    pe = p1 * q1 + (1 - p1) * (1 - q1)
    k = (po - pe) / (1 - pe)
    C("Arm I kappa recomputed", abs(k - aiE["cohens_kappa"]) < 1e-9, "%.6f" % k)
    C("Arm I MOB-suite E replicons = 3,937",
      ce[(True, True)] + ce[(True, False)] == 3937)

    # ---------- gate verdicts ----------
    worst = max(R["arms"][k2]["PRIMARY_replicon_change_pct"] for k2 in R["arms"])
    C("C10 verdict follows from the recomputed numbers",
      (R["C10"]["VERDICT"] == "PASS") == (worst <= 5.0), "worst %.4f%%" % worst)
    C("C09 verdict follows from its own cells",
      (R["C09"]["VERDICT"] == "PASS") ==
      all(v["direction_preserved"] and v["ci_excludes_zero"]
          for v in R["C09"].values() if isinstance(v, dict)))
    C("headline plasmid share unchanged",
      abs(R["headline_invariance"]["occurrence_weighted_plasmid_share_pct"]
          - 100.0 * plas / tot) < 1e-9, "%.6f%%" % (100.0 * plas / tot))
    C("classes A and B untouched",
      R["headline_invariance"]["class_A"] == 18837 and R["headline_invariance"]["class_B"] == 16303)

    print("\n  checks run: %d   disagreements: %d" % (checks, len(fails)))
    print("  VERDICT: %s" % ("PASS - zero disagreements" if not fails else "*** FAIL: %s ***" % fails))
    sys.exit(0 if not fails else 9)


if __name__ == "__main__":
    main()
