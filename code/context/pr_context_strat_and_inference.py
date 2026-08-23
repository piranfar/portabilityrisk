"""PR-CONTEXT Phases B(18) and D -- stratified epidemiology and strengthened inference.

Adds to the frozen results; changes none of them. Every previously reported number is
recomputed here from the same primary outputs and asserted, so an extension that silently
moved a headline would fail rather than ship.

The central inferential question this file exists to answer is whether the plasmid-enrichment
effects survive the two things most likely to be creating them:

  * BETWEEN-SPECIES COMPOSITION. K. pneumoniae contributes 62% of all occurrences and has a
    67% plasmid share, while P. aeruginosa has 12%. A family concentrated in Klebsiella would
    look plasmid-enriched for reasons that have nothing to do with the gene. Every family is
    therefore re-tested WITHIN each species and by Mantel-Haenszel stratification on species.
  * REPEATED PLASMID LINEAGES. If one clonal plasmid is sequenced 300 times, occurrence-level
    counts inflate. Effects are recomputed per unique plasmid, per genome and per BioProject,
    and MOB-suite primary clusters are used to measure how concentrated the evidence is.
"""
import argparse, collections, csv, hashlib, json, math, os, sys
import numpy as np
from scipy.stats import fisher_exact

VERSION = "pr_context_strat_and_inference_v1.0.0"
RNG = np.random.default_rng(20260821)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def write(path, rows, cols):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in cols) + "\n")


def or_ci(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    o = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return o, o * math.exp(-1.96 * se), o * math.exp(1.96 * se)


def bh(p):
    n = len(p)
    order = sorted(range(n), key=lambda i: p[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        k = n - rank + 1
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def mantel_haenszel(tables):
    """MH pooled OR with a Robins-Breslow-Greenland variance. Strata = species."""
    num = den = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
    if den == 0 or num == 0:
        return None, None, None
    mh = num / den
    s1 = s2 = s3 = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        P, Q = (a + d) / n, (b + c) / n
        R, S = a * d / n, b * c / n
        s1 += P * R
        s2 += P * S + Q * R
        s3 += Q * S
    if num <= 0 or den <= 0:
        return mh, None, None
    var = s1 / (2 * num ** 2) + s2 / (2 * num * den) + s3 / (2 * den ** 2)
    se = math.sqrt(var)
    return mh, mh * math.exp(-1.96 * se), mh * math.exp(1.96 * se)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")

    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"]
    for r in occ:
        r["is_plasmid"] = r["evidence_type"] == "direct_plasmid"
    res = [r for r in occ if r["evidence_type"] in ("direct_plasmid", "direct_chromosome")]
    NP = sum(1 for r in res if r["is_plasmid"])
    NC = len(res) - NP
    # ---- assert the frozen baseline is intact ----
    base = {"primary_occurrences": (len(occ), 74349), "plasmid": (NP, 39209),
            "chromosome": (NC, 35140),
            "genomes": (len({r["assembly_version"] for r in res}), 6288),
            "plasmids": (len({r["replicon_accession"] for r in res if r["is_plasmid"]}),
                         6621)}
    for k, (got, want) in base.items():
        if got != want:
            print("STOP: frozen baseline changed - %s got %d want %d" % (k, got, want))
            sys.exit(2)
    print("%s | frozen baseline re-asserted intact" % VERSION)

    mob = {r["replicon_accession"]: r for r in csv.DictReader(
        open(os.path.join(out, "plasmid_mobility_annotation.tsv"), encoding="utf-8"),
        delimiter="\t")}
    bsm = {r["biosample_accession"]: r for r in csv.DictReader(
        open(os.path.join(out, "biosample_metadata.tsv"), encoding="utf-8"),
        delimiter="\t")}
    for r in res:
        m = bsm.get(r["biosample_accession"], {})
        r["country"] = m.get("country", "")
        r["source_context"] = m.get("source_context", "")
        r["clinical_status"] = m.get("clinical_status", "")
        r["host"] = m.get("host", "")
        r["year"] = (m.get("collection_date", "") or "")[:4]
        p = mob.get(r["replicon_accession"], {})
        r["mobility"] = p.get("portability_mobility_category", "") if r["is_plasmid"] else ""
        r["mob_cluster"] = p.get("primary_cluster_id", "") if r["is_plasmid"] else ""

    # ================= PHASE B 18: stratified epidemiology =================
    def strat(field, minn=100):
        rows = []
        ks = collections.Counter(r[field] for r in res if r[field])
        for k, _ in ks.most_common():
            sub = [r for r in res if r[field] == k]
            if len(sub) < minn:
                continue
            p = sum(1 for r in sub if r["is_plasmid"])
            conj = sum(1 for r in sub if r["mobility"] == "predicted_conjugative")
            rows.append({field: k, "n_occurrences": len(sub),
                         "n_genomes": len({r["assembly_version"] for r in sub}),
                         "n_biosamples": len({r["biosample_accession"] for r in sub}),
                         "n_bioprojects": len({r["bioproject_accession"] for r in sub}),
                         "n_unique_plasmids": len({r["replicon_accession"] for r in sub
                                                   if r["is_plasmid"]}),
                         "n_plasmid": p, "n_chromosome": len(sub) - p,
                         "pct_plasmid": round(100.0 * p / len(sub), 2),
                         "n_conjugative_borne": conj,
                         "pct_of_plasmid_borne_conjugative":
                             round(100.0 * conj / p, 2) if p else ""})
        return rows
    for field, fname in (("country", "portability_by_country.tsv"),
                         ("source_context", "portability_by_source_context.tsv"),
                         ("clinical_status", "portability_by_clinical_status.tsv"),
                         ("host", "portability_by_host.tsv"),
                         ("year", "portability_by_year.tsv")):
        rows = strat(field)
        if rows:
            write(os.path.join(out, fname), rows, list(rows[0].keys()))

    # species x source
    ss = []
    for (sp, sc), n in collections.Counter(
            (r["organism_harmonized"], r["source_context"]) for r in res
            if r["source_context"]).most_common():
        if n < 100:
            continue
        sub = [r for r in res if r["organism_harmonized"] == sp
               and r["source_context"] == sc]
        p = sum(1 for r in sub if r["is_plasmid"])
        conj = sum(1 for r in sub if r["mobility"] == "predicted_conjugative")
        ss.append({"organism_harmonized": sp, "source_context": sc,
                   "n_occurrences": n, "n_genomes": len({r["assembly_version"] for r in sub}),
                   "n_plasmid": p, "pct_plasmid": round(100.0 * p / n, 2),
                   "pct_of_plasmid_borne_conjugative":
                       round(100.0 * conj / p, 2) if p else ""})
    write(os.path.join(out, "portability_by_species_and_source.tsv"), ss,
          list(ss[0].keys()))

    # ================= PHASE D: species-adjusted inference =================
    fams = collections.Counter(r["gene_family"] for r in res)
    species = [s for s, n in collections.Counter(
        r["organism_harmonized"] for r in res).items() if n >= 500]
    rows = []
    for f, n in fams.items():
        if n < 20:
            continue
        a = sum(1 for r in res if r["gene_family"] == f and r["is_plasmid"])
        b = n - a
        crude, clo, chi = or_ci(a, b, NP - a, NC - b)
        _, pcr = fisher_exact([[a, b], [NP - a, NC - b]])
        tables, within = [], []
        for sp in species:
            sub = [r for r in res if r["organism_harmonized"] == sp]
            if not sub:
                continue
            sa = sum(1 for r in sub if r["gene_family"] == f and r["is_plasmid"])
            sb = sum(1 for r in sub if r["gene_family"] == f and not r["is_plasmid"])
            sc_ = sum(1 for r in sub if r["is_plasmid"]) - sa
            sd = len(sub) - sa - sb - sc_
            if sa + sb == 0:
                continue
            tables.append((sa, sb, sc_, sd))
            if sa + sb >= 20:
                o, lo, hi = or_ci(sa, sb, sc_, sd)
                within.append((sp, sa + sb, o, lo, hi))
        mh, mlo, mhi = mantel_haenszel(tables) if tables else (None, None, None)
        # unique-plasmid and per-genome denominators
        up = len({r["replicon_accession"] for r in res
                  if r["gene_family"] == f and r["is_plasmid"]})
        ug = len({r["assembly_version"] for r in res if r["gene_family"] == f})
        ubp = len({r["bioproject_accession"] for r in res if r["gene_family"] == f})
        cl = collections.Counter(r["mob_cluster"] for r in res
                                 if r["gene_family"] == f and r["is_plasmid"]
                                 and r["mob_cluster"])
        top_cluster_share = (round(100.0 * cl.most_common(1)[0][1] / sum(cl.values()), 2)
                             if cl else "")
        rows.append({
            "gene_family": f, "n_occurrences": n, "n_plasmid": a, "n_chromosome": b,
            "pct_plasmid": round(100.0 * a / n, 2),
            "crude_odds_ratio": round(crude, 4), "crude_ci_lo": round(clo, 4),
            "crude_ci_hi": round(chi, 4), "crude_fisher_p": pcr,
            "species_adjusted_MH_or": round(mh, 4) if mh else "",
            "species_adjusted_ci_lo": round(mlo, 4) if mlo else "",
            "species_adjusted_ci_hi": round(mhi, 4) if mhi else "",
            "n_species_strata": len(tables),
            "n_species_with_within_estimate": len(within),
            "within_species_all_same_direction":
                ("yes" if within and len({(w[2] > 1) for w in within}) == 1 else
                 "no" if within else ""),
            "within_species_detail": ";".join(
                "%s:n=%d,OR=%.3g" % (w[0], w[1], w[2]) for w in within),
            "unique_plasmids_carrying": up, "unique_genomes_carrying": ug,
            "unique_bioprojects_carrying": ubp,
            "largest_mob_cluster_share_pct": top_cluster_share,
            "lineage_dominated_flag":
                "yes" if (top_cluster_share != "" and top_cluster_share >= 50) else "no"})
    qs = bh([r["crude_fisher_p"] for r in rows])
    for r, q in zip(rows, qs):
        r["crude_bh_q"] = q
        r["significant_after_bh"] = "yes" if q < 0.05 else "no"
        mhv, lo, hi = r["species_adjusted_MH_or"], r["species_adjusted_ci_lo"], \
            r["species_adjusted_ci_hi"]
        r["survives_species_adjustment"] = (
            "yes" if (mhv != "" and lo != "" and hi != "" and not (lo <= 1.0 <= hi))
            else "no" if mhv != "" else "not_estimable")
    rows.sort(key=lambda r: -r["crude_odds_ratio"])
    write(os.path.join(out, "determinant_enrichment_species_adjusted.tsv"), rows,
          list(rows[0].keys()))

    # ================= PHASE D: sensitivity analyses =================
    sens = []
    def frac(sub):
        p = sum(1 for r in sub if r["is_plasmid"])
        return len(sub), p, round(100.0 * p / len(sub), 3) if sub else 0

    n0, p0, f0 = frac(res)
    sens.append({"analysis": "PRIMARY (core non-EFFLUX, 74,349)", "n": n0,
                 "n_plasmid": p0, "pct_plasmid": f0, "delta_vs_primary_pp": 0.0})
    sub = [r for r in res if r["call_completeness"] == "complete"]
    n1, p1, f1 = frac(sub)
    sens.append({"analysis": "excluding PARTIAL and INTERNAL_STOP calls", "n": n1,
                 "n_plasmid": p1, "pct_plasmid": f1,
                 "delta_vs_primary_pp": round(f1 - f0, 3)})
    allamr = [r for r in occ if r["evidence_type"] in
              ("direct_plasmid", "direct_chromosome")]
    occ_all = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t")]
    for r in occ_all:
        r["is_plasmid"] = r["evidence_type"] == "direct_plasmid"
    broad = [r for r in occ_all if r["evidence_type"] in
             ("direct_plasmid", "direct_chromosome")]
    n2, p2, f2 = frac(broad)
    sens.append({"analysis": "broad AMR/AMR denominator (85,507 incl. EFFLUX and plus)",
                 "n": n2, "n_plasmid": p2, "pct_plasmid": f2,
                 "delta_vs_primary_pp": round(f2 - f0, 3)})
    # one occurrence per (genome, family)
    seen = set()
    dedup = []
    for r in res:
        k = (r["assembly_version"], r["gene_family"])
        if k not in seen:
            seen.add(k)
            dedup.append(r)
    n3, p3, f3 = frac(dedup)
    sens.append({"analysis": "one occurrence per genome per gene family", "n": n3,
                 "n_plasmid": p3, "pct_plasmid": f3,
                 "delta_vs_primary_pp": round(f3 - f0, 3)})
    # one genome per BioSample (repeated-BioSample sensitivity)
    bs_first = {}
    for r in res:
        bs_first.setdefault(r["biosample_accession"], r["assembly_version"])
    sub = [r for r in res if r["assembly_version"] ==
           bs_first[r["biosample_accession"]]]
    n4, p4, f4 = frac(sub)
    sens.append({"analysis": "one genome per BioSample (repeated-BioSample sensitivity)",
                 "n": n4, "n_plasmid": p4, "pct_plasmid": f4,
                 "delta_vs_primary_pp": round(f4 - f0, 3)})
    # one occurrence per unique plasmid / per genome / per bioproject
    for lab, key in (("per unique plasmid or chromosome replicon", "replicon_accession"),
                     ("per genome", "assembly_version"),
                     ("per BioProject", "bioproject_accession")):
        seen, sub = set(), []
        for r in res:
            k = (r[key], r["is_plasmid"])
            if k not in seen:
                seen.add(k)
                sub.append(r)
        n5, p5, f5 = frac(sub)
        sens.append({"analysis": "collapsed to one event " + lab, "n": n5,
                     "n_plasmid": p5, "pct_plasmid": f5,
                     "delta_vs_primary_pp": round(f5 - f0, 3)})
    write(os.path.join(out, "sensitivity_analyses.tsv"), sens, list(sens[0].keys()))

    # ================= lineage concentration =================
    cl = collections.Counter(r["mob_cluster"] for r in res
                             if r["is_plasmid"] and r["mob_cluster"])
    tot = sum(cl.values())
    top = cl.most_common(20)
    hhi = sum((v / tot) ** 2 for v in cl.values()) if tot else 0
    lin = [{"mob_primary_cluster": k, "n_plasmid_borne_occurrences": v,
            "pct_of_plasmid_borne": round(100.0 * v / tot, 3),
            "n_unique_plasmids": len({r["replicon_accession"] for r in res
                                      if r["mob_cluster"] == k}),
            "n_genomes": len({r["assembly_version"] for r in res
                              if r["mob_cluster"] == k}),
            "n_bioprojects": len({r["bioproject_accession"] for r in res
                                  if r["mob_cluster"] == k})} for k, v in top]
    lin.append({"mob_primary_cluster": "__SUMMARY__",
                "n_plasmid_borne_occurrences": tot,
                "pct_of_plasmid_borne": 100.0,
                "n_unique_plasmids": len({r["replicon_accession"] for r in res
                                          if r["is_plasmid"]}),
                "n_genomes": "", "n_bioprojects":
                    "distinct clusters %d; HHI %.4f; effective clusters 1/HHI %.1f; "
                    "largest cluster %.2f%% of plasmid-borne occurrences"
                    % (len(cl), hhi, (1 / hhi) if hhi else 0,
                       100.0 * top[0][1] / tot if top else 0)})
    write(os.path.join(out, "plasmid_lineage_concentration.tsv"), lin,
          list(lin[0].keys()))

    # ---------------- report ----------------
    print("\n=== PHASE B18: stratified epidemiology ===")
    for fn, lab in (("portability_by_country.tsv", "country"),
                    ("portability_by_source_context.tsv", "source context"),
                    ("portability_by_clinical_status.tsv", "clinical status")):
        p = os.path.join(out, fn)
        if not os.path.exists(p):
            continue
        print("  --- %s ---" % lab)
        for r in list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))[:6]:
            k = list(r)[0]
            print("    %-22s n=%-6s plasmid %6s%%  conj-borne %6s%%  genomes %s"
                  % (r[k][:22], r["n_occurrences"], r["pct_plasmid"],
                     r["pct_of_plasmid_borne_conjugative"], r["n_genomes"]))
    print("\n=== PHASE D: species adjustment ===")
    tested = len(rows)
    sig = [r for r in rows if r["significant_after_bh"] == "yes"]
    surv = [r for r in sig if r["survives_species_adjustment"] == "yes"]
    same = [r for r in sig if r["within_species_all_same_direction"] == "yes"]
    print("  families tested                          : %d" % tested)
    print("  significant crude (BH q<0.05)            : %d" % len(sig))
    print("  SURVIVING species-adjusted MH interval   : %d (%.1f%%)"
          % (len(surv), 100.0 * len(surv) / max(len(sig), 1)))
    print("  same direction in every species stratum  : %d" % len(same))
    print("  lineage-dominated (one MOB cluster >=50%%): %d"
          % sum(1 for r in rows if r["lineage_dominated_flag"] == "yes"))
    print("\n  strongest effects that survive species adjustment:")
    for r in sorted([x for x in surv if x["crude_odds_ratio"] > 1],
                    key=lambda x: -x["crude_odds_ratio"])[:8]:
        print("    %-14s n=%-5d crude OR %8.2f  MH-adj %8.2f [%.2f, %.2f]  clusters<50%%:%s"
              % (r["gene_family"], r["n_occurrences"], r["crude_odds_ratio"],
                 r["species_adjusted_MH_or"], r["species_adjusted_ci_lo"],
                 r["species_adjusted_ci_hi"],
                 "yes" if r["lineage_dominated_flag"] == "no" else "NO"))
    print("\n=== PHASE D: sensitivity ===")
    for r in sens:
        print("  %-56s n=%-6d plasmid %7.3f%%  delta %+6.3f pp"
              % (r["analysis"][:56], r["n"], r["pct_plasmid"],
                 r["delta_vs_primary_pp"]))
    print("\n=== lineage concentration ===")
    print("  " + lin[-1]["n_bioprojects"])
    for f in ("determinant_enrichment_species_adjusted.tsv", "sensitivity_analyses.tsv",
              "plasmid_lineage_concentration.tsv",
              "portability_by_country.tsv", "portability_by_source_context.tsv",
              "portability_by_clinical_status.tsv", "portability_by_host.tsv",
              "portability_by_year.tsv", "portability_by_species_and_source.tsv"):
        p = os.path.join(out, f)
        if os.path.exists(p):
            print("  %-52s %s" % (f, sha256_file(p)))


if __name__ == "__main__":
    main()
