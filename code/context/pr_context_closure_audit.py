"""PR-CONTEXT closure audit -- read only.

Builds a per-file manifest carrying relative_path, bytes, sha256, artefact_class and
required_for_reproduction, so that "is everything here?" can be answered per artefact rather
than as a single aggregate. Runs identically on the laptop and on the server; the two are
compared afterwards.

Nothing is written outside the audit directory, no analysis is rerun, and no frozen result is
touched.
"""
import argparse, csv, hashlib, json, os, re, sys

VERSION = "pr_context_closure_audit_v1.0.0"

# artefact_class, required_for_reproduction
RULES = [
 (re.compile(r"^out/FROZEN_.*\.json$"),                ("frozen_protocol", "yes")),
 (re.compile(r"^out/frozen_portability_class_.*json$"), ("frozen_protocol", "yes")),
 (re.compile(r"^out/future_model_join_schema\.json$"), ("integration_schema", "yes")),
 (re.compile(r"^out/evidence_layer_dictionary\.tsv$"), ("integration_schema", "yes")),
 (re.compile(r"^out/independent_verification.*"),      ("verification", "yes")),
 (re.compile(r"^out/headline_recomputation\.tsv$"),    ("verification", "yes")),
 (re.compile(r"^out/portability_class_reconciliation\.tsv$"), ("verification", "yes")),
 (re.compile(r"^out/.*receipts?\.(tsv|json)$"),        ("receipt", "yes")),
 (re.compile(r"^out/DELIVERABLES_MANIFEST\.tsv$"),     ("manifest", "yes")),
 (re.compile(r"^out/environment.*"),                   ("environment", "yes")),
 (re.compile(r"^out/commands\.txt$"),                  ("commands", "yes")),
 (re.compile(r"^out/.*\.(tsv|json|txt)$"),             ("deliverable", "yes")),
 (re.compile(r"^scripts/.*\.py$"),                     ("script", "yes")),
 (re.compile(r"^report/PORTABILITYRISK_.*"),           ("final_deliverable", "yes")),
 (re.compile(r"^report/figures/.*\.png$"),             ("figure", "yes")),
 (re.compile(r"^closure/manifest_.*\.tsv$"),           ("manifest", "yes")),
 (re.compile(r"^closure/REPORT_FACTS\.json$"),         ("manifest", "yes")),
 (re.compile(r"^closure/BUNDLE_.*"),                   ("manifest", "yes")),
 (re.compile(r"^closure/.*\.tar\.gz$"),               ("bundle", "no")),
 (re.compile(r"^logs/.*"),                             ("log", "no")),
 (re.compile(r"^plasmid_fasta/.*\.fna$"),
  ("sequence_cache_plasmid", "no - reconstructable from plasmid_retrieval_receipts.tsv "
                             "(accession.version + raw_sha256)")),
 (re.compile(r"^window_fasta/.*\.fna$"),
  ("sequence_cache_window", "no - reconstructable from shared_context_blocks.tsv "
                            "(replicon + coordinates + extracted_sequence_sha256)")),
 (re.compile(r"^seqreport_cache/.*"),   ("metadata_cache", "no - re-fetchable from NCBI")),
 (re.compile(r"^biosample_cache/.*"),   ("metadata_cache", "no - re-fetchable from NCBI")),
 (re.compile(r"^topo_cache/.*"),        ("metadata_cache", "no - re-fetchable from NCBI")),
 (re.compile(r"^mob_out/.*"),           ("tool_output_raw", "no - regenerable by MOB-suite")),
 (re.compile(r"^mge_work/.*"),          ("tool_work_scratch", "no")),
 (re.compile(r"^mobdb/.*"),             ("tool_database", "no - re-installable via mob_init")),
]
SKIP_DIRS = {"mge_work", "mobdb", "ref_pkg", "pkg"}


def classify(rel):
    for pat, val in RULES:
        if pat.match(rel):
            return val
    return ("other", "no")


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--scripts", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    a = ap.parse_args()

    rows = []
    for base, sub in [(a.root, None)] + ([(a.scripts, "scripts")] if a.scripts else []):
        if not os.path.isdir(base):
            continue
        for r, dirs, files in os.walk(base):
            dirs[:] = [d for d in dirs
                       if os.path.relpath(os.path.join(r, d), base).split(os.sep)[0]
                       not in SKIP_DIRS]
            for f in files:
                p = os.path.join(r, f)
                rel = os.path.relpath(p, base).replace("\\", "/")
                if sub:
                    rel = sub + "/" + rel
                # A manifest cannot consistently describe itself: the server writes its
                # copy, the laptop then receives it, and the two enumerations can never
                # agree on that one file. Same class as the bundle receipt. Excluded.
                if (rel.startswith("audit_manifest") or rel.endswith(".tmp")
                        or rel in ("closure/manifest_local.tsv",
                                   "closure/manifest_server.tsv")):
                    continue
                cls, req = classify(rel)
                try:
                    rows.append({"relative_path": rel, "bytes": os.path.getsize(p),
                                 "sha256": sha256_file(p), "artefact_class": cls,
                                 "required_for_reproduction": req})
                except OSError:
                    continue
    rows.sort(key=lambda r: r["relative_path"])
    cols = ["relative_path", "bytes", "sha256", "artefact_class",
            "required_for_reproduction"]
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    import collections
    cc = collections.Counter(r["artefact_class"] for r in rows)
    bb = collections.Counter()
    for r in rows:
        bb[r["artefact_class"]] += r["bytes"]
    req = [r for r in rows if r["required_for_reproduction"] == "yes"]
    agg = hashlib.sha256(("\n".join("%s\t%d\t%s" % (r["relative_path"], r["bytes"],
                                                    r["sha256"]) for r in rows)
                          + "\n").encode()).hexdigest()
    aggreq = hashlib.sha256(("\n".join("%s\t%d\t%s" % (r["relative_path"], r["bytes"],
                                                       r["sha256"]) for r in req)
                             + "\n").encode()).hexdigest()
    print("%s [%s]" % (VERSION, a.label))
    print("  files %d | bytes %d (%.3f GB)"
          % (len(rows), sum(r["bytes"] for r in rows),
             sum(r["bytes"] for r in rows) / 1e9))
    print("  required-for-reproduction files %d | bytes %d (%.1f MB)"
          % (len(req), sum(r["bytes"] for r in req),
             sum(r["bytes"] for r in req) / 2 ** 20))
    for k in sorted(cc, key=lambda x: -bb[x]):
        print("    %-26s %6d files  %9.3f MB" % (k, cc[k], bb[k] / 2 ** 20))
    print("  AGGREGATE(all)      %s" % agg)
    print("  AGGREGATE(required) %s" % aggreq)


if __name__ == "__main__":
    main()
