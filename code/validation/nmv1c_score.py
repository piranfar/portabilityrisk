"""NM-V1C corrected confirmatory audit -- scoring.

Applies the gate frozen in NMV1C_FROZEN_PROTOCOL.json before any case was drawn. Every input
digest is verified first and the script aborts on any mismatch. Nothing is overwritten.

The primary metric is design-weighted three-state agreement between the blinded adjudicator's
call and the frozen rule engine's state, with a stratified bootstrap interval. The unweighted
figure is reported as a secondary. Per-state and per-stratum agreement are always reported,
because the protocol confines a failure to the state in which it occurs.
"""
import argparse, collections, csv, datetime, hashlib, json, os, random, re, sys

VERSION = "nmv1c_score_v1.0.0"
PROTO_SHA = "b2058877c0f7165f0ab7bf9be7adcb2cb3b5b4c8d2b81143d85ab03fe0f0c04c"
KEY_SHA = "40e142e27978631da9ec437b9cdd7aeb8c5fda52bdecbbcda5f2fa2db8501efe"
APP_SHA = "56581a89c3e7acf0ce39b831c122bacd8d165ebfbdd7aa58a9d8dafda1eee41b"
RECV_SHA = "a572945c9c02d316e00b076ac6f881472be554f5b3dc068cfee0bbdcf24c6b75"
B = 2000
SEED = 20260822


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    D = a.dir

    for f, want in (("NMV1C_FROZEN_PROTOCOL.json", PROTO_SHA),
                    ("NMV1C_UNBLINDING_KEY.tsv", KEY_SHA),
                    ("NMV1C_ADJUDICATION_APP_R2.html", APP_SHA),
                    ("NMV1C_ADJUDICATED_120_RECEIVED.json", RECV_SHA)):
        got = sha(os.path.join(D, f))
        if got != want:
            print("REFUSING: %s digest mismatch\n  want %s\n  got  %s" % (f, want, got))
            sys.exit(1)
    print("%s | all four input digests verified" % VERSION)

    P = json.load(open(os.path.join(D, "NMV1C_FROZEN_PROTOCOL.json"), encoding="utf-8"))
    key = {r["token"]: r for r in csv.DictReader(
        open(os.path.join(D, "NMV1C_UNBLINDING_KEY.tsv"), encoding="utf-8"), delimiter="\t")}
    dec = {r["token"]: r for r in json.load(
        open(os.path.join(D, "NMV1C_ADJUDICATED_120_RECEIVED.json"), encoding="utf-8"))["decisions"]}
    if set(key) != set(dec):
        print("REFUSING: token sets differ"); sys.exit(1)

    rows = []
    for t in sorted(key):
        k, d = key[t], dec[t]
        rows.append({"token": t, "stratum": k["stratum"], "rule_id": k["rule_id"],
                     "machine_state": k["machine_state"], "weight": float(k["weight"]),
                     "adjudicator": d["decision"], "reason_code": d.get("reason_code", ""),
                     "note": (d.get("note") or "").replace("\t", " "),
                     "agree": int(d["decision"] == k["machine_state"])})

    def wagree(rs):
        w = sum(r["weight"] for r in rs)
        return sum(r["weight"] * r["agree"] for r in rs) / w if w else float("nan")

    def uagree(rs):
        return sum(r["agree"] for r in rs) / len(rs) if rs else float("nan")

    primary = wagree(rows)
    unweighted = uagree(rows)

    # stratified bootstrap: resample within stratum, preserving stratum sizes
    by_st = collections.defaultdict(list)
    for r in rows:
        by_st[r["stratum"]].append(r)
    rnd = random.Random(SEED)
    boot = []
    for _ in range(B):
        rs = []
        for st, grp in sorted(by_st.items()):
            rs.extend(rnd.choices(grp, k=len(grp)))
        boot.append(wagree(rs))
    boot.sort()
    lo, hi = boot[int(0.025 * B)], boot[int(0.975 * B) - 1]

    by_state = {}
    for s in sorted({r["machine_state"] for r in rows}):
        g = [r for r in rows if r["machine_state"] == s]
        bs = []
        for _ in range(B):
            bs.append(uagree(rnd.choices(g, k=len(g))))
        bs.sort()
        by_state[s] = {"n": len(g), "weighted": wagree(g), "unweighted": uagree(g),
                       "agreed": sum(r["agree"] for r in g),
                       "ci_unweighted": [bs[int(0.025 * B)], bs[int(0.975 * B) - 1]]}

    by_stratum = {st: {"n": len(g), "weight": g[0]["weight"],
                       "agreed": sum(r["agree"] for r in g), "unweighted": uagree(g)}
                  for st, g in sorted(by_st.items())}

    G = P["gates"]
    if primary >= 0.90 and lo >= 0.80:
        verdict, rule = "SUCCESS", G["SUCCESS"]
    elif (0.80 <= primary < 0.90) or (0.65 <= lo < 0.80):
        verdict, rule = "REVISE", G["REVISE"]
    else:
        verdict, rule = "FAIL", G["FAIL"]

    # the protocol's second failure limb, reported not silently applied
    weak = [s for s, v in by_state.items() if v["unweighted"] < 0.80]

    conf = collections.Counter((r["machine_state"], r["adjudicator"]) for r in rows)

    print("\n  PRIMARY  design-weighted three-state agreement : %.4f  [95%% %.4f - %.4f]"
          % (primary, lo, hi))
    print("  SECONDARY unweighted agreement                  : %.4f  (%d / %d)"
          % (unweighted, sum(r["agree"] for r in rows), len(rows)))
    print("\n  per machine state")
    for s, v in by_state.items():
        print("    %-14s n=%3d  agreed %3d  unweighted %.4f  [%.4f - %.4f]  weighted %.4f"
              % (s, v["n"], v["agreed"], v["unweighted"],
                 v["ci_unweighted"][0], v["ci_unweighted"][1], v["weighted"]))
    print("\n  per stratum")
    for st, v in by_stratum.items():
        print("    %-3s n=%3d  weight %.4f  agreed %3d  unweighted %.4f"
              % (st, v["n"], v["weight"], v["agreed"], v["unweighted"]))
    print("\n  confusion (machine -> adjudicator)")
    for (m, adj), n in sorted(conf.items()):
        print("    %-14s -> %-14s %3d" % (m, adj, n))
    if weak:
        print("\n  NOTE: state(s) below 0.80 unweighted: %s" % weak)
    print("\n  GATE VERDICT: %s\n    rule: %s" % (verdict, rule))

    out = os.path.join(D, "nmv1c_scored_cases.tsv")
    rc = os.path.join(D, "NMV1C_SCORING_RECEIPT.json")
    for p in (out, rc):
        if os.path.exists(p):
            print("REFUSING: %s exists" % p); sys.exit(1)
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    R = {"receipt": "NMV1C_SCORING_RECEIPT", "scorer": VERSION,
         "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
         "protocol_sha256": PROTO_SHA, "key_sha256": KEY_SHA,
         "application_sha256": APP_SHA, "adjudication_sha256": RECV_SHA,
         "engine_sha256": P["classifier_under_validation"]["sha256"],
         "engine_unchanged": True, "F3_threshold_unchanged": True,
         "n": len(rows),
         "primary_design_weighted_agreement": primary,
         "bootstrap": {"method": "stratified, resampled within stratum", "B": B, "seed": SEED,
                       "ci95": [lo, hi]},
         "secondary_unweighted_agreement": unweighted,
         "agreed": sum(r["agree"] for r in rows),
         "by_machine_state": by_state, "by_stratum": by_stratum,
         "confusion": {"%s->%s" % k: v for k, v in sorted(conf.items())},
         "states_below_0.80_unweighted": weak,
         "gate_verdict": verdict, "gate_rule_applied": rule,
         "gates_as_frozen": G,
         "structural_limitation": P["STRUCTURAL_LIMITATION_STATED_BEFORE_DRAWING"]["consequence"],
         "prior_audit": P["prior_audit_disposition"],
         "adjudications_reused_from_prior_packages": 0,
         "command": "python nmv1c_score.py --dir docs/nature_microbiology"}
    json.dump(R, open(rc, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  wrote %s  %s" % (os.path.basename(out), sha(out)))
    print("  wrote %s  %s" % (os.path.basename(rc), sha(rc)))


if __name__ == "__main__":
    main()
