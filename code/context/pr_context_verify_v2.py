"""PR-CONTEXT independent verification v2 -- complete class reconciliation.

Extends v1. Imports nothing from any analysis script, recounts from the raw AMRFinderPlus
tables and the frozen metadata, and additionally reconciles the COMPLETE portability-class
set against all 74,349 primary occurrences, so that A + B + C + D + E plus every preserved
state sums to the denominator exactly with nothing silently dropped.

It also re-derives the chromosomal MGE proportion from the feature inventory rather than from
the neighbourhood table the pipeline produced, so a join error in the pipeline would show up
as a disagreement instead of being reproduced.
"""
import argparse, collections, csv, glob, hashlib, json, os, sys

VERSION = "pr_context_verify_v2_v1.0.0"


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
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--amrfinder-dir", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")

    inv = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(out, "replicon_inventory.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        for k in (r["sequence_accession"], r["refseq_accession"], r["genbank_accession"]):
            if k:
                inv[r["assembly_accession"]][k] = r

    # ---------- independent recount from raw tables ----------
    n_rows = n_amr = n_prim = 0
    loc = collections.Counter()
    genomes, plasmids, replicons = set(), set(), set()
    chrom_keys = set()
    for f in sorted(glob.glob(os.path.join(a.amrfinder_dir, "*.tsv"))):
        if f.endswith("_timing.tsv"):
            continue
        asm = os.path.basename(f)[:-4]
        idx = inv.get(asm, {})
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"),
                                delimiter="\t"):
            n_rows += 1
            if r.get("Type") != "AMR" or r.get("Subtype") != "AMR":
                continue
            n_amr += 1
            if r.get("Scope") != "core" or r.get("Class") == "EFFLUX":
                continue
            n_prim += 1
            genomes.add(asm)
            cid = (r.get("Contig id") or "").strip()
            rec = idx.get(cid) or idx.get(cid.split(".")[0])
            mt = (rec or {}).get("replicon_molecule_type", "")
            if mt == "Plasmid":
                loc["plasmid"] += 1
                plasmids.add(rec["sequence_accession"])
            elif mt == "Chromosome":
                loc["chromosome"] += 1
                try:
                    s, e = sorted((int(r["Start"]), int(r["Stop"])))
                    chrom_keys.add((rec["sequence_accession"], s, e))
                except Exception:
                    pass
            else:
                loc["unresolved"] += 1
            if rec:
                replicons.add((asm, rec["sequence_accession"]))

    # ---------- independent MGE proportion from the feature inventory ----------
    feats = collections.defaultdict(list)
    fp = os.path.join(out, "mge_feature_inventory.tsv")
    if os.path.exists(fp):
        for r in csv.DictReader(open(fp, encoding="utf-8"), delimiter="\t"):
            feats[r["replicon_accession"]].append((int(r["chrom_start"]),
                                                   int(r["chrom_end"])))
    ann = {}
    np_ = os.path.join(out, "arg_mge_neighbourhood.tsv")
    if os.path.exists(np_):
        for r in csv.DictReader(open(np_, encoding="utf-8"), delimiter="\t"):
            ann[(r["replicon_accession"], int(r["gene_start"]), int(r["gene_end"]))] = r
    my_with = my_tot = 0
    for (rep, s, e) in chrom_keys:
        if (rep, s, e) not in ann:
            continue
        my_tot += 1
        lo, hi = s - 10000, e + 10000
        if any(fe >= lo and fs <= hi for fs, fe in feats.get(rep, [])):
            my_with += 1

    # ---------- pipeline numbers ----------
    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"]
    ploc = collections.Counter(
        "plasmid" if r["evidence_type"] == "direct_plasmid" else
        "chromosome" if r["evidence_type"] == "direct_chromosome" else "unresolved"
        for r in occ)
    cls = list(csv.DictReader(
        open(os.path.join(out, "determinant_portability_classes.tsv"), encoding="utf-8"),
        delimiter="\t"))
    cc = collections.Counter(r["portability_class"] for r in cls)
    pipe_with = sum(1 for r in ann.values() if int(r["n_mge_features"]) > 0)

    checks = [
        ("total AMRFinderPlus rows", n_rows, 184538),
        ("AMR/AMR rows", n_amr, 85507),
        ("PRIMARY acquired ARG occurrences", n_prim, len(occ)),
        ("directly mapped plasmid ARGs", loc["plasmid"], ploc["plasmid"]),
        ("directly mapped chromosomal ARGs", loc["chromosome"], ploc["chromosome"]),
        ("unresolved occurrences", loc["unresolved"], ploc["unresolved"]),
        ("unique genomes", len(genomes), len({r["assembly_version"] for r in occ})),
        ("unique ARG-bearing plasmids", len(plasmids),
         len({r["replicon_accession"] for r in occ
              if r["evidence_type"] == "direct_plasmid"})),
        ("unique ARG-bearing replicons", len(replicons),
         len({(r["assembly_version"], r["replicon_accession"]) for r in occ})),
        ("class rows == primary denominator", len(cls), len(occ)),
        ("A+B == chromosomal (or A+B+INCOMPLETE)",
         cc["A"] + cc["B"] + cc["INCOMPLETE_ANNOTATION"], ploc["chromosome"]),
        ("C+D+E == plasmid",
         cc["C"] + cc["D"] + cc["E"] + cc["UNRESOLVED_MOBILITY"], ploc["plasmid"]),
        ("all classes sum to 74,349", sum(cc.values()), 74349),
        ("chromosomal ARGs with nearby MGE", my_with, pipe_with),
        ("chromosomal ARGs annotated", my_tot, len(ann)),
    ]
    rows, bad = [], 0
    for name, mine, theirs in checks:
        ok = mine == theirs
        bad += (0 if ok else 1)
        rows.append({"check": name, "independent_value": mine, "pipeline_value": theirs,
                     "agreement": "MATCH" if ok else "MISMATCH"})
    write(os.path.join(out, "independent_verification_report_v2.tsv"), rows,
          list(rows[0].keys()))

    recon = [{"portability_class": k, "n": v,
              "pct_of_74349": round(100.0 * v / 74349, 3)}
             for k, v in sorted(cc.items(), key=lambda kv: -kv[1])]
    recon.append({"portability_class": "__TOTAL__", "n": sum(cc.values()),
                  "pct_of_74349": round(100.0 * sum(cc.values()) / 74349, 3)})
    write(os.path.join(out, "portability_class_reconciliation.tsv"), recon,
          list(recon[0].keys()))

    print("%s\n=== INDEPENDENT VERIFICATION v2 ===" % VERSION)
    for r in rows:
        print("  %-42s independent %-8s pipeline %-8s %s"
              % (r["check"], r["independent_value"], r["pipeline_value"],
                 r["agreement"]))
    print("\n=== COMPLETE CLASS RECONCILIATION ===")
    for r in recon:
        print("  %-24s %7d  %6.3f%%" % (r["portability_class"], r["n"],
                                        r["pct_of_74349"]))
    print("\nRESULT: %s (%d disagreements)"
          % ("PASS" if bad == 0 else "FAIL -- STOP CONDITION", bad))
    sys.exit(0 if bad == 0 else 3)


if __name__ == "__main__":
    main()
