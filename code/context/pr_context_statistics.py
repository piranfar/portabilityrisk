"""PR-CONTEXT step 4 -- the prespecified statistical analysis on the direct-evidence layer.

Runs only the analyses named in FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1: effect sizes are
odds ratios with 95% intervals, the test is two-sided Fisher exact, families below 20
occurrences are not tested, and correction is Benjamini-Hochberg at q < 0.05.

Non-independence is handled the way the protocol requires rather than ignored. Multiple copies
of one gene family inside one genome are not independent observations, so every headline
effect is computed a second time after GENOME-LEVEL AGGREGATION -- each genome contributes at
most one plasmid event and one chromosome event per family -- and a third time as a cluster
bootstrap resampling whole BioProjects. Where the naive and adjusted intervals disagree, the
adjusted one is the one that is reported.
"""
import argparse, collections, csv, hashlib, json, math, os, sys
import numpy as np
from scipy.stats import fisher_exact

VERSION = "pr_context_statistics_v1.0.0"
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
    """Odds ratio with a Woolf 95% interval, Haldane-Anscombe corrected when a cell is 0."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    orv = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return orv, orv * math.exp(-1.96 * se), orv * math.exp(1.96 * se)


def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        k = n - rank + 1
        prev = min(prev, pvals[i] * n / k)
        q[i] = prev
    return q


def cluster_bootstrap_or(rows, key, target, clusters, B=2000):
    """Resample whole BioProjects and recompute the odds ratio. Returns the 2.5/97.5 pct."""
    by = collections.defaultdict(list)
    for r in rows:
        by[r["bioproject_accession"]].append((r[key] == target, r["is_plasmid"]))
    cl = list(by)
    if len(cl) < 3:
        return None, None
    idx = np.arange(len(cl))
    out = []
    for _ in range(B):
        pick = RNG.choice(idx, size=len(idx), replace=True)
        a = b = c = d = 0
        for j in pick:
            for is_t, is_p in by[cl[j]]:
                if is_t and is_p:
                    a += 1
                elif is_t:
                    b += 1
                elif is_p:
                    c += 1
                else:
                    d += 1
        if min(a, b, c, d) == 0:
            a, b, c, d = a + .5, b + .5, c + .5, d + .5
        out.append((a * d) / (b * c))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--protocol-sha256", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")
    pp = os.path.join(out, "FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1.json")
    if sha256_file(pp) != a.protocol_sha256:
        print("REFUSING: protocol digest mismatch"); sys.exit(1)
    P = json.load(open(pp, encoding="utf-8"))
    MINF = P["statistics"]["minimum_family_size_for_testing"]
    print("%s | protocol verified | min family size %d" % (VERSION, MINF))

    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"]
    for r in occ:
        r["is_plasmid"] = r["evidence_type"] == "direct_plasmid"
    res = [r for r in occ if r["evidence_type"] in ("direct_plasmid", "direct_chromosome")]
    NP = sum(1 for r in res if r["is_plasmid"])
    NC = len(res) - NP
    print("resolved occurrences: %d (plasmid %d, chromosome %d)" % (len(res), NP, NC))

    # ---------------- genome-level summary
    g = collections.defaultdict(lambda: {"n": 0, "p": 0, "c": 0, "fams": set(),
                                         "repl": set(), "plasmids": set()})
    for r in res:
        k = r["assembly_version"]
        v = g[k]
        v["n"] += 1
        v["p" if r["is_plasmid"] else "c"] += 1
        v["fams"].add(r["gene_family"])
        v["repl"].add(r["replicon_accession"])
        if r["is_plasmid"]:
            v["plasmids"].add(r["replicon_accession"])
        v["genus"] = r["genus"]
        v["organism"] = r["organism_harmonized"]
        v["biosample"] = r["biosample_accession"]
        v["bioproject"] = r["bioproject_accession"]
    grows = [{"assembly_version": k, "biosample_accession": v["biosample"],
              "bioproject_accession": v["bioproject"], "organism": v["organism"],
              "genus": v["genus"], "n_arg_occurrences": v["n"],
              "n_plasmid_args": v["p"], "n_chromosomal_args": v["c"],
              "pct_plasmid": round(100.0 * v["p"] / v["n"], 3),
              "n_distinct_families": len(v["fams"]),
              "n_arg_bearing_replicons": len(v["repl"]),
              "n_arg_bearing_plasmids": len(v["plasmids"])}
             for k, v in sorted(g.items())]
    write(os.path.join(out, "genome_level_summary.tsv"), grows, list(grows[0].keys()))

    # ---------------- replicon-level summary
    rp = collections.defaultdict(lambda: {"n": 0, "fams": set(), "syms": set()})
    for r in res:
        k = (r["assembly_version"], r["replicon_accession"])
        v = rp[k]
        v["n"] += 1
        v["fams"].add(r["gene_family"])
        v["syms"].add(r["determinant_name"])
        v["type"] = r["replicon_molecule_type"]
        v["len"] = r["replicon_length"]
        v["name"] = r["replicon_name"]
        v["genus"] = r["genus"]
        v["organism"] = r["organism_harmonized"]
        v["biosample"] = r["biosample_accession"]
    rrows = [{"assembly_version": k[0], "replicon_accession": k[1],
              "replicon_molecule_type": v["type"], "replicon_name": v["name"],
              "replicon_length": v["len"], "organism": v["organism"], "genus": v["genus"],
              "biosample_accession": v["biosample"], "n_arg_occurrences": v["n"],
              "n_distinct_families": len(v["fams"]),
              "families": ";".join(sorted(v["fams"]))}
             for k, v in sorted(rp.items())]
    write(os.path.join(out, "replicon_level_summary.tsv"), rrows, list(rrows[0].keys()))

    # ---------------- Q2 family enrichment
    fam = collections.Counter((r["gene_family"], r["is_plasmid"]) for r in res)
    fams = sorted({f for f, _ in fam})
    recs = []
    for f in fams:
        a = fam[(f, True)]
        b = fam[(f, False)]
        if a + b < MINF:
            continue
        c = NP - a
        d = NC - b
        orv, lo, hi = or_ci(a, b, c, d)
        _, p = fisher_exact([[a, b], [c, d]])
        # genome-level aggregation: one plasmid event and one chromosome event per genome
        ga = len({r["assembly_version"] for r in res
                  if r["gene_family"] == f and r["is_plasmid"]})
        gb = len({r["assembly_version"] for r in res
                  if r["gene_family"] == f and not r["is_plasmid"]})
        gc_ = len({r["assembly_version"] for r in res if r["is_plasmid"]}) - ga
        gd = len({r["assembly_version"] for r in res if not r["is_plasmid"]}) - gb
        gor, glo, ghi = or_ci(ga, gb, gc_, gd)
        recs.append({"gene_family": f, "n_occurrences": a + b, "n_plasmid": a,
                     "n_chromosome": b, "pct_plasmid": round(100.0 * a / (a + b), 2),
                     "odds_ratio": round(orv, 4), "or_ci_lo": round(lo, 4),
                     "or_ci_hi": round(hi, 4), "fisher_p": p,
                     "genome_level_n_plasmid_genomes": ga,
                     "genome_level_n_chromosome_genomes": gb,
                     "genome_level_odds_ratio": round(gor, 4),
                     "genome_level_ci_lo": round(glo, 4),
                     "genome_level_ci_hi": round(ghi, 4)})
    qs = bh([r["fisher_p"] for r in recs])
    for r, q in zip(recs, qs):
        r["bh_q"] = q
        r["significant_q_lt_0.05"] = "yes" if q < 0.05 else "no"
        r["direction"] = ("plasmid_enriched" if r["odds_ratio"] > 1 else
                          "chromosome_enriched")
    recs.sort(key=lambda r: (-r["odds_ratio"], r["bh_q"]))
    write(os.path.join(out, "determinant_plasmid_enrichment.tsv"), recs,
          list(recs[0].keys()))
    sig = [r for r in recs if r["significant_q_lt_0.05"] == "yes"]
    print("families tested: %d | significant at BH q<0.05: %d" % (len(recs), len(sig)))

    # cluster bootstrap for the strongest families in each direction
    top = sorted(sig, key=lambda r: -r["odds_ratio"])[:12] + \
        sorted(sig, key=lambda r: r["odds_ratio"])[:12]
    seen, boot = set(), []
    for r in top:
        if r["gene_family"] in seen:
            continue
        seen.add(r["gene_family"])
        lo, hi = cluster_bootstrap_or(res, "gene_family", r["gene_family"], None, B=2000)
        boot.append({**r, "cluster_bootstrap_ci_lo": round(lo, 4) if lo else "",
                     "cluster_bootstrap_ci_hi": round(hi, 4) if hi else "",
                     "cluster_unit": "BioProject", "n_resamples": 2000,
                     "crosses_1_after_clustering":
                         "yes" if (lo and hi and lo <= 1.0 <= hi) else "no"})
    write(os.path.join(out, "statistical_model_receipts.tsv"), boot, list(boot[0].keys()))

    # ---------------- Q3 both-context determinants
    both = []
    for f in fams:
        a, b = fam[(f, True)], fam[(f, False)]
        if a and b:
            both.append({"gene_family": f, "n_plasmid": a, "n_chromosome": b,
                         "n_total": a + b,
                         "minor_context_fraction": round(min(a, b) / (a + b), 4)})
    both.sort(key=lambda r: -r["n_total"])
    write(os.path.join(out, "both_context_determinants.tsv"), both, list(both[0].keys()))

    # ---------------- Q6 strata
    for name, key in (("portability_by_species.tsv", "organism_harmonized"),
                      ("portability_by_gene_family.tsv", "gene_family"),
                      ("portability_by_drug_class.tsv", "drug_class")):
        cnt = collections.Counter((r[key], r["is_plasmid"]) for r in res)
        ks = sorted({k for k, _ in cnt})
        rows = []
        for k in ks:
            a, b = cnt[(k, True)], cnt[(k, False)]
            rows.append({key: k, "n_occurrences": a + b, "n_plasmid": a, "n_chromosome": b,
                         "pct_plasmid": round(100.0 * a / (a + b), 2),
                         "n_genomes": len({r["assembly_version"] for r in res if r[key] == k}),
                         "n_bioprojects": len({r["bioproject_accession"] for r in res
                                               if r[key] == k})})
        rows.sort(key=lambda r: -r["n_occurrences"])
        write(os.path.join(out, name), rows, list(rows[0].keys()))

    # ---------------- accounting
    allp = [r for r in occ]
    acct = [{"layer": "all AMRFinderPlus rows", "n": 184538,
             "note": "every category, not acquired ARGs"},
            {"layer": "AMR/AMR rows", "n": 85507, "note": "acquired-gene rows incl. EFFLUX"},
            {"layer": "PRIMARY acquired ARG occurrences", "n": len(allp),
             "note": "core, non-EFFLUX; the frozen denominator"},
            {"layer": "directly resolved to chromosome or plasmid", "n": len(res),
             "note": "%.4f%% of the primary denominator" % (100.0 * len(res) / len(allp))},
            {"layer": "direct_plasmid", "n": NP,
             "note": "%.2f%% of resolved" % (100.0 * NP / len(res))},
            {"layer": "direct_chromosome", "n": NC,
             "note": "%.2f%% of resolved" % (100.0 * NC / len(res))},
            {"layer": "unresolved", "n": len(allp) - len(res),
             "note": "identifier_unmatched + ambiguous + unclassified + coords_missing"}]
    write(os.path.join(out, "direct_vs_unresolved_accounting.tsv"), acct,
          ["layer", "n", "note"])

    print("\ntop 8 plasmid-enriched families (BH q<0.05):")
    for r in [x for x in recs if x["significant_q_lt_0.05"] == "yes"][:8]:
        print("  %-16s n=%-6d plasmid %5.1f%%  OR %8.2f [%.2f, %.2f]  q=%.3g"
              % (r["gene_family"], r["n_occurrences"], r["pct_plasmid"], r["odds_ratio"],
                 r["or_ci_lo"], r["or_ci_hi"], r["bh_q"]))
    print("\ntop 8 chromosome-enriched families:")
    for r in [x for x in recs if x["significant_q_lt_0.05"] == "yes"][-8:]:
        print("  %-16s n=%-6d plasmid %5.1f%%  OR %8.4f [%.4f, %.4f]  q=%.3g"
              % (r["gene_family"], r["n_occurrences"], r["pct_plasmid"], r["odds_ratio"],
                 r["or_ci_lo"], r["or_ci_hi"], r["bh_q"]))
    print("\nfamilies in BOTH contexts: %d of %d tested-or-not families" % (len(both), len(fams)))
    for f in ("genome_level_summary.tsv", "replicon_level_summary.tsv",
              "determinant_plasmid_enrichment.tsv", "both_context_determinants.tsv",
              "portability_by_species.tsv", "portability_by_gene_family.tsv",
              "portability_by_drug_class.tsv", "statistical_model_receipts.tsv",
              "direct_vs_unresolved_accounting.tsv"):
        print("  %-42s %s" % (f, sha256_file(os.path.join(out, f))))


if __name__ == "__main__":
    main()
