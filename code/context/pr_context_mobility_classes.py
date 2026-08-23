"""PR-CONTEXT step 7 -- plasmid mobility, and the frozen portability-evidence classes.

The class definitions are written to disk and hashed BEFORE any class is assigned.

Two wording rules from the task are load-bearing and are implemented literally rather than
paraphrased:

  * a plasmid is never called non-mobilizable because a marker database returned no match.
    MOB-suite's "non-mobilizable" is recorded as
    nonconjugative_or_no_mobility_markers_detected, which is a statement about detection.
  * a plasmid is never called conjugative from a single isolated transposase. Conjugative
    here requires MOB-suite to report a relaxase AND a mate-pair-formation type, which is a
    multi-component requirement by construction.

The classes are an EVIDENCE-BASED PORTABILITY CLASS. They are not a calibrated risk score and
not a transmission probability: no transmission outcome was measured anywhere in this study.
"""
import argparse, collections, csv, datetime, glob, hashlib, json, os, sys

VERSION = "pr_context_mobility_classes_v1.0.0"

DEFS = {
 "framework": "evidence-based portability class",
 "version": "1.0.0",
 "not_a_risk_score":
   "These classes describe WHAT WAS OBSERVED in sequence, nothing else. They are not a "
   "calibrated clinical risk score and must never be described as a measured transmission "
   "probability: no transmission outcome exists anywhere in this dataset.",
 "dimensions_are_separate": [
   "direct genomic location -- NCBI replicon designation on a closed replicon",
   "nearby MGE signature -- sequence annotation in a fixed window",
   "plasmid mobilization evidence -- MOB-suite relaxase / oriT",
   "conjugation-associated evidence -- MOB-suite relaxase plus mate-pair formation",
   "data completeness -- recorded explicitly, never imputed"],
 "classes": {
   "A": {"definition": "direct chromosomal location AND no prespecified MGE feature "
                       "detected in the frozen window",
         "requires": ["evidence_type == direct_chromosome",
                      "MGE annotation COMPLETED for this window",
                      "zero qualifying MGE features"]},
   "B": {"definition": "direct chromosomal location AND at least one prespecified MGE "
                       "feature in the frozen window",
         "requires": ["evidence_type == direct_chromosome",
                      "MGE annotation COMPLETED", ">= 1 qualifying MGE feature"]},
   "C": {"definition": "direct plasmid location, insufficient genomic evidence for "
                       "mobilization",
         "requires": ["evidence_type == direct_plasmid",
                      "MOB-suite typed the replicon",
                      "no relaxase AND no oriT detected"]},
   "D": {"definition": "direct plasmid location on a replicon with evidence consistent "
                       "with mobilization",
         "requires": ["evidence_type == direct_plasmid",
                      "relaxase OR oriT detected",
                      "NOT the full conjugation combination required for class E"]},
   "E": {"definition": "direct plasmid location on a replicon with combined evidence "
                       "consistent with conjugation",
         "requires": ["evidence_type == direct_plasmid",
                      "relaxase detected AND mate-pair-formation type detected"]}},
 "preserved_states": {
   "UNRESOLVED_LOCATION": "no direct chromosome/plasmid designation for the occurrence",
   "UNRESOLVED_MOBILITY": "plasmid location established but MOB-suite produced no usable "
                          "record for that replicon",
   "INCOMPLETE_ANNOTATION": "the dimension needed to separate two classes has not been "
                            "computed yet for this occurrence. Never merged into A or B."},
 "mobility_categories": {
   "predicted_conjugative": "relaxase AND mpf_type reported by MOB-suite",
   "predicted_mobilizable": "relaxase OR oriT reported, without the full conjugation set",
   "nonconjugative_or_no_mobility_markers_detected":
     "MOB-suite reported neither relaxase nor oriT. A DETECTION statement about a marker "
     "database, NOT a claim that the plasmid cannot be mobilized.",
   "mobility_unresolved": "no usable MOB-suite record for the replicon"},
 "evidence_type_vocabulary": {
   "detected_sequence_feature": "a sequence match was found at the stated thresholds",
   "predicted_functional_mobility": "a database-based prediction about capability",
   "experimentally_demonstrated_mobility":
     "NEVER assigned in this study. No transfer experiment exists here."},
 "software": {"mob_typer_version": "3.1.9",
              "mob_typer_sha256":
                "43408d898356cecc872589f33e8f8a2b551012f6f46eb4edb8bf190d7a66f954",
              "database_recursive_sha256":
                "aae85a3bc4a51b1cc859d83e9b1c8b1abd337c9b5ee1b29d5c94af619c6624d5",
              "database_size": "1.6 GB", "blast": "2.15.0", "mash": "2.3",
              "invocation": "mob_typer --multi -n 1 -i <chunk.fna> -o <out> -d <db>",
              "thresholds": "MOB-suite 3.1.9 defaults, unmodified",
              "result_nature": "database-based prediction, not direct experimental evidence"},
}


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


def present(v):
    return bool(v) and v.strip() not in ("", "-", "NA", "nan", "none", "None")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")
    dp = os.path.join(out, "frozen_portability_class_definitions.json")
    if not os.path.exists(dp):
        DEFS["created_utc"] = datetime.datetime.now(
            datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        DEFS["created_before_any_class_was_assigned"] = True
        json.dump(DEFS, open(dp, "w", encoding="utf-8", newline="\n"), indent=2)
    print("%s\nclass definitions frozen: %s" % (VERSION, sha256_file(dp)))

    # ---------------- MOB-suite aggregation
    mob = {}
    for f in sorted(glob.glob(os.path.join(a.root, "mob_out", "*.txt"))):
        for r in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t"):
            sid = (r.get("sample_id") or "").split()[0]
            mob[sid] = r
    print("MOB-suite records: %d" % len(mob))

    prows = []
    for acc, r in sorted(mob.items()):
        rel = present(r.get("relaxase_type(s)"))
        mpf = present(r.get("mpf_type"))
        ori = present(r.get("orit_type(s)"))
        rep = present(r.get("rep_type(s)"))
        if rel and mpf:
            cat = "predicted_conjugative"
        elif rel or ori:
            cat = "predicted_mobilizable"
        else:
            cat = "nonconjugative_or_no_mobility_markers_detected"
        prows.append({
            "replicon_accession": acc, "plasmid_length": r.get("size", ""),
            "gc": r.get("gc", ""), "md5": r.get("md5", ""),
            "rep_types": r.get("rep_type(s)", ""),
            "rep_type_accessions": r.get("rep_type_accession(s)", ""),
            "relaxase_types": r.get("relaxase_type(s)", ""),
            "mpf_type": r.get("mpf_type", ""),
            "orit_types": r.get("orit_type(s)", ""),
            "mob_predicted_mobility": r.get("predicted_mobility", ""),
            "portability_mobility_category": cat,
            "has_relaxase": "yes" if rel else "no",
            "has_mpf": "yes" if mpf else "no",
            "has_orit": "yes" if ori else "no",
            "has_replicon_type": "yes" if rep else "no",
            "primary_cluster_id": r.get("primary_cluster_id", ""),
            "predicted_host_range_rank": r.get("predicted_host_range_overall_rank", ""),
            "predicted_host_range_name": r.get("predicted_host_range_overall_name", ""),
            "evidence_type": "predicted_functional_mobility",
            "annotation_source": "MOB-suite 3.1.9 mob_typer, defaults"})
    write(os.path.join(out, "plasmid_mobility_annotation.tsv"), prows,
          list(prows[0].keys()))
    write(os.path.join(out, "plasmid_replicon_typing.tsv"),
          [{k: p[k] for k in ("replicon_accession", "plasmid_length", "rep_types",
                              "rep_type_accessions", "has_replicon_type",
                              "primary_cluster_id", "predicted_host_range_rank",
                              "predicted_host_range_name")} for p in prows],
          ["replicon_accession", "plasmid_length", "rep_types", "rep_type_accessions",
           "has_replicon_type", "primary_cluster_id", "predicted_host_range_rank",
           "predicted_host_range_name"])

    # ---------------- occurrences -> classes
    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"]
    mge_done = os.path.exists(os.path.join(out, "arg_mge_neighbourhood.tsv"))
    mge = {}
    if mge_done:
        for r in csv.DictReader(open(os.path.join(out, "arg_mge_neighbourhood.tsv"),
                                     encoding="utf-8"), delimiter="\t"):
            mge[(r["assembly_version"], r["replicon_accession"], r["gene_start"])] = r
    print("MGE annotation available: %s" % mge_done)

    crows, comp = [], []
    for r in occ:
        et = r["evidence_type"]
        pl = mob.get(r["replicon_accession"])
        cls = mob_cat = ""
        if et == "direct_plasmid":
            if pl is None:
                cls, mob_cat = "UNRESOLVED_MOBILITY", "mobility_unresolved"
            else:
                rel = present(pl.get("relaxase_type(s)"))
                mpf = present(pl.get("mpf_type"))
                ori = present(pl.get("orit_type(s)"))
                if rel and mpf:
                    cls, mob_cat = "E", "predicted_conjugative"
                elif rel or ori:
                    cls, mob_cat = "D", "predicted_mobilizable"
                else:
                    cls, mob_cat = ("C",
                                    "nonconjugative_or_no_mobility_markers_detected")
        elif et == "direct_chromosome":
            k = (r["assembly_version"], r["replicon_accession"], str(r["gene_start"]))
            if not mge_done or k not in mge:
                cls, mob_cat = "INCOMPLETE_ANNOTATION", "not_applicable_chromosomal"
            else:
                n = int(mge[k].get("n_mge_features", 0) or 0)
                cls = "B" if n > 0 else "A"
                mob_cat = "not_applicable_chromosomal"
        else:
            cls, mob_cat = "UNRESOLVED_LOCATION", "location_unresolved"
        crows.append({
            "assembly_version": r["assembly_version"],
            "biosample_accession": r["biosample_accession"],
            "bioproject_accession": r["bioproject_accession"],
            "organism_harmonized": r["organism_harmonized"], "genus": r["genus"],
            "replicon_accession": r["replicon_accession"],
            "replicon_molecule_type": r["replicon_molecule_type"],
            "determinant_name": r["determinant_name"], "gene_family": r["gene_family"],
            "drug_class": r["drug_class"], "gene_start": r["gene_start"],
            "gene_end": r["gene_end"], "strand": r["strand"],
            "location_evidence": et,
            "portability_class": cls,
            "mobility_category": mob_cat,
            "evidence_layer_location": "direct_closed_replicon",
            "evidence_layer_mobility": ("predicted_plasmid_mobility"
                                        if et == "direct_plasmid" else ""),
            "evidence_layer_mge": ("mge_sequence_annotation" if mge_done
                                   and et == "direct_chromosome" else ""),
            "plasmidcall_predicted_location": ""})
        comp.append({"assembly_version": r["assembly_version"],
                     "replicon_accession": r["replicon_accession"],
                     "determinant_name": r["determinant_name"],
                     "gene_start": r["gene_start"],
                     "component_location": et,
                     "component_relaxase": (pl or {}).get("relaxase_type(s)", ""),
                     "component_mpf": (pl or {}).get("mpf_type", ""),
                     "component_orit": (pl or {}).get("orit_type(s)", ""),
                     "component_rep_type": (pl or {}).get("rep_type(s)", ""),
                     "component_mge": "PENDING" if not mge_done else "",
                     "portability_class": cls})
    write(os.path.join(out, "determinant_portability_classes.tsv"), crows,
          list(crows[0].keys()))
    write(os.path.join(out, "portability_class_evidence_components.tsv"), comp,
          list(comp[0].keys()))

    cc = collections.Counter(r["portability_class"] for r in crows)
    print("\n=== PORTABILITY CLASS DISTRIBUTION over %d acquired ARG occurrences ==="
          % len(crows))
    for k in ("A", "B", "C", "D", "E", "UNRESOLVED_LOCATION", "UNRESOLVED_MOBILITY",
              "INCOMPLETE_ANNOTATION"):
        n = cc.get(k, 0)
        if n:
            print("  %-24s %7d  %6.2f%%" % (k, n, 100.0 * n / len(crows)))
    pl_tot = sum(cc.get(k, 0) for k in ("C", "D", "E")) + cc.get("UNRESOLVED_MOBILITY", 0)
    if pl_tot:
        print("\n  --- Q4: of %d plasmid-borne ARG occurrences:" % pl_tot)
        for k, lab in (("E", "conjugation-consistent"), ("D", "mobilization-consistent"),
                       ("C", "no mobility markers detected"),
                       ("UNRESOLVED_MOBILITY", "unresolved")):
            n = cc.get(k, 0)
            print("      %-32s %7d  %6.2f%%" % (lab, n, 100.0 * n / pl_tot))
    mc = collections.Counter(p["portability_mobility_category"] for p in prows)
    print("\n=== unique ARG-bearing plasmids: %d ===" % len(prows))
    for k, v in mc.most_common():
        print("  %-50s %5d  %5.1f%%" % (k, v, 100.0 * v / len(prows)))
    rep = collections.Counter()
    for p in prows:
        for t in (p["rep_types"] or "").split(","):
            if t.strip():
                rep[t.strip()] += 1
    print("\ntop replicon/Inc types:")
    for k, v in rep.most_common(10):
        print("  %-24s %5d" % (k, v))
    miss = [{"metric": "plasmids typed", "n": len(prows)},
            {"metric": "plasmids with a replicon type", "n": sum(1 for p in prows if p["has_replicon_type"] == "yes")},
            {"metric": "plasmids with a relaxase", "n": sum(1 for p in prows if p["has_relaxase"] == "yes")},
            {"metric": "plasmids with an MPF type", "n": sum(1 for p in prows if p["has_mpf"] == "yes")},
            {"metric": "plasmids with an oriT", "n": sum(1 for p in prows if p["has_orit"] == "yes")},
            {"metric": "ARG occurrences awaiting MGE annotation",
             "n": cc.get("INCOMPLETE_ANNOTATION", 0)}]
    write(os.path.join(out, "plasmid_mobility_missingness.tsv"), miss, ["metric", "n"])
    for f in ("frozen_portability_class_definitions.json",
              "plasmid_mobility_annotation.tsv", "plasmid_replicon_typing.tsv",
              "determinant_portability_classes.tsv",
              "portability_class_evidence_components.tsv",
              "plasmid_mobility_missingness.tsv"):
        print("  %-48s %s" % (f, sha256_file(os.path.join(out, f))))


if __name__ == "__main__":
    main()
