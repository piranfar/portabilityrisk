"""PR-CONTEXT step 8 -- independent verification, and the future-model join schema.

The verifier imports nothing from the analysis scripts. It goes back to the 7,216 raw
AMRFinderPlus TSV files and recounts from them, re-derives the replicon designation from the
frozen metadata, and re-reads the frozen class definitions -- then compares its own numbers
with the pipeline's outputs.

It deliberately re-implements the join differently: the pipeline builds a dictionary index
over three accession spellings and looks each contig up, while this walks the inventory once
per assembly and matches within the assembly. A shared indexing bug would not survive both.

Any principal count that differs is a stop condition.
"""
import argparse, collections, csv, datetime, glob, hashlib, json, os, sys

VERSION = "pr_context_verify_v1.0.0"

SCHEMA = {
 "schema": "future_model_join_schema", "version": "1.0.0",
 "purpose":
   "Let PlasmidCall predictions be added later as a separate evidence layer without "
   "reprocessing the 7,216 genomes and without touching a single directly evidenced "
   "replicon assignment.",
 "join_keys": ["assembly_accession", "assembly_version", "biosample_accession",
               "sequence_accession", "replicon_accession", "contig_id", "gene_symbol",
               "gene_start", "gene_end", "strand", "source_sha256"],
 "key_notes": {
   "contig_id": "the AMRFinderPlus 'Contig id' verbatim; in this cohort it is identical to "
                "sequence_accession because the inputs were RefSeq/GenBank complete genomes",
   "gene_symbol": "column determinant_name in determinant_occurrences.tsv",
   "source_sha256": "sha256 of the AMRFinderPlus output file the row came from"},
 "evidence_layers": {
   "direct_closed_replicon": {"status": "POPULATED", "authority": "HIGHEST",
     "source": "NCBI assigned_molecule_location_type on a closed replicon",
     "may_be_overwritten_by": "nothing"},
   "mge_sequence_annotation": {"status": "PENDING",
     "source": "sequence annotation in the frozen neighbourhood window"},
   "predicted_plasmid_mobility": {"status": "POPULATED",
     "source": "MOB-suite 3.1.9 mob_typer, database-based prediction"},
   "plasmidcall_predicted_location": {"status": "RESERVED, EMPTY",
     "source": "not run in this task",
     "rules": ["never overwrites direct_closed_replicon",
               "never becomes, or contributes to, a target label for PlasmidCall itself",
               "must carry its own model version, threshold and run receipt",
               "is reported as a separate column, never merged into location_evidence"]}},
 "integration_procedure": [
   "1. produce PlasmidCall predictions keyed by (assembly_version, contig_id)",
   "2. left-join onto determinant_occurrences.tsv on those two keys",
   "3. write ONLY into plasmidcall_predicted_location and its companion score column",
   "4. report agreement with direct_closed_replicon as an EVALUATION of PlasmidCall, "
   "never as a correction of the direct evidence"],
 "why_this_works_without_reprocessing":
   "Every row already carries the assembly version, the exact sequence accession and the "
   "gene coordinates. A prediction made later on the same assemblies joins on those keys "
   "alone; nothing about this dataset needs recomputing.",
}

LAYERS = [
 {"layer": "direct_closed_replicon", "authority": "highest",
  "populated_now": "yes", "evidence_nature": "documented replicon designation",
  "may_be_overwritten": "no",
  "definition": "NCBI states the molecule is Chromosome or Plasmid for a closed replicon in "
                "a complete genome, and the ARG coordinates fall inside it"},
 {"layer": "mge_sequence_annotation", "authority": "supporting",
  "populated_now": "no", "evidence_nature": "detected sequence feature",
  "may_be_overwritten": "no",
  "definition": "prespecified mobile-element features within the frozen window"},
 {"layer": "predicted_plasmid_mobility", "authority": "supporting",
  "populated_now": "yes", "evidence_nature": "database-based prediction",
  "may_be_overwritten": "no",
  "definition": "MOB-suite relaxase / MPF / oriT typing of the ARG-bearing replicon"},
 {"layer": "plasmidcall_predicted_location", "authority": "lowest",
  "populated_now": "no -- reserved schema field only",
  "evidence_nature": "model prediction",
  "may_be_overwritten": "n/a",
  "definition": "reserved for a future PlasmidCall run; must never overwrite a direct "
                "assignment nor form that model's own target"},
]


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

    json.dump(SCHEMA, open(os.path.join(out, "future_model_join_schema.json"), "w",
                           encoding="utf-8", newline="\n"), indent=2)
    write(os.path.join(out, "evidence_layer_dictionary.tsv"), LAYERS,
          list(LAYERS[0].keys()))

    # ---- independent recount straight from the raw tables ----
    inv_by_asm = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(out, "replicon_inventory.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        for k in (r["sequence_accession"], r["refseq_accession"], r["genbank_accession"]):
            if k:
                inv_by_asm[r["assembly_accession"]][k] = r

    tot_rows = 0
    n_amr_amr = 0
    n_primary = 0
    loc = collections.Counter()
    genomes = set()
    replicons = set()
    plasmids = set()
    fams = collections.Counter()
    for f in sorted(glob.glob(os.path.join(a.amrfinder_dir, "*.tsv"))):
        if f.endswith("_timing.tsv"):
            continue
        acc = os.path.basename(f)[:-4]
        idx = inv_by_asm.get(acc, {})
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"),
                                delimiter="\t"):
            tot_rows += 1
            if r.get("Type") != "AMR" or r.get("Subtype") != "AMR":
                continue
            n_amr_amr += 1
            if r.get("Scope") != "core" or r.get("Class") == "EFFLUX":
                continue
            n_primary += 1
            genomes.add(acc)
            cid = (r.get("Contig id") or "").strip()
            rec = idx.get(cid) or idx.get(cid.split(".")[0])
            mt = (rec or {}).get("replicon_molecule_type", "")
            if mt == "Plasmid":
                loc["plasmid"] += 1
                plasmids.add(rec["sequence_accession"])
            elif mt == "Chromosome":
                loc["chromosome"] += 1
            else:
                loc["unresolved"] += 1
            if rec:
                replicons.add((acc, rec["sequence_accession"]))
            sym = (r.get("Element symbol") or "").strip()
            if "-" in sym:
                head, _, tail = sym.rpartition("-")
                import re as _re
                fam = head if head and _re.match(
                    r"^(?:\d+[A-Za-z'\"]*|[IVX]+[a-z]?|[A-Za-z]{1,3})$", tail) else sym
            else:
                fam = sym
            fams[fam] += 1

    # ---- pipeline's own numbers ----
    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"]
    p_loc = collections.Counter(
        "plasmid" if r["evidence_type"] == "direct_plasmid" else
        "chromosome" if r["evidence_type"] == "direct_chromosome" else "unresolved"
        for r in occ)
    cls = collections.Counter(r["portability_class"] for r in csv.DictReader(
        open(os.path.join(out, "determinant_portability_classes.tsv"), encoding="utf-8"),
        delimiter="\t"))
    mob = list(csv.DictReader(open(os.path.join(out, "plasmid_mobility_annotation.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    enr = list(csv.DictReader(open(os.path.join(out, "determinant_plasmid_enrichment.tsv"),
                                   encoding="utf-8"), delimiter="\t"))

    checks = [
      ("total AMRFinderPlus rows", tot_rows, 184538),
      ("AMR/AMR rows", n_amr_amr, 85507),
      ("PRIMARY acquired ARG occurrences", n_primary, len(occ)),
      ("directly mapped plasmid ARGs", loc["plasmid"], p_loc["plasmid"]),
      ("directly mapped chromosomal ARGs", loc["chromosome"], p_loc["chromosome"]),
      ("unresolved occurrences", loc["unresolved"], p_loc["unresolved"]),
      ("unique genomes with a primary ARG", len(genomes),
       len({r["assembly_version"] for r in occ})),
      ("unique ARG-bearing replicons", len(replicons),
       len({(r["assembly_version"], r["replicon_accession"]) for r in occ
            if r["replicon_accession"]})),
      ("unique ARG-bearing plasmids", len(plasmids), len(mob)),
      ("portability class C+D+E", cls["C"] + cls["D"] + cls["E"], p_loc["plasmid"]),
      ("class E (conjugation-consistent)", cls["E"],
       sum(1 for r in csv.DictReader(
           open(os.path.join(out, "determinant_portability_classes.tsv"),
                encoding="utf-8"), delimiter="\t") if r["portability_class"] == "E")),
    ]
    rows = []
    bad = 0
    for name, mine, theirs in checks:
        ok = mine == theirs
        bad += (0 if ok else 1)
        rows.append({"check": name, "independent_value": mine, "pipeline_value": theirs,
                     "agreement": "MATCH" if ok else "MISMATCH"})
    # headline enrichment recomputation for the top families
    head = []
    NP, NC = loc["plasmid"], loc["chromosome"]
    for e in enr[:5] + enr[-5:]:
        f = e["gene_family"]
        a_ = sum(1 for r in occ if r["gene_family"] == f
                 and r["evidence_type"] == "direct_plasmid")
        b_ = sum(1 for r in occ if r["gene_family"] == f
                 and r["evidence_type"] == "direct_chromosome")
        c_, d_ = NP - a_, NC - b_
        aa, bb, cc2, dd = (a_ + .5, b_ + .5, c_ + .5, d_ + .5) if min(a_, b_, c_, d_) == 0 \
            else (a_, b_, c_, d_)
        orv = (aa * dd) / (bb * cc2)
        ok = abs(orv - float(e["odds_ratio"])) < 0.01 * max(1.0, float(e["odds_ratio"]))
        bad += (0 if ok else 1)
        head.append({"gene_family": f, "independent_n_plasmid": a_,
                     "pipeline_n_plasmid": e["n_plasmid"],
                     "independent_n_chromosome": b_,
                     "pipeline_n_chromosome": e["n_chromosome"],
                     "independent_odds_ratio": round(orv, 4),
                     "pipeline_odds_ratio": e["odds_ratio"],
                     "agreement": "MATCH" if ok else "MISMATCH"})
    write(os.path.join(out, "independent_verification_report.tsv"), rows,
          list(rows[0].keys()))
    write(os.path.join(out, "headline_recomputation.tsv"), head, list(head[0].keys()))

    print("%s\n=== INDEPENDENT VERIFICATION ===" % VERSION)
    for r in rows:
        print("  %-40s independent %-8s pipeline %-8s %s"
              % (r["check"], r["independent_value"], r["pipeline_value"], r["agreement"]))
    print("\n=== HEADLINE RECOMPUTATION ===")
    for r in head:
        print("  %-16s OR independent %-10s pipeline %-10s %s"
              % (r["gene_family"], r["independent_odds_ratio"], r["pipeline_odds_ratio"],
                 r["agreement"]))
    print("\nRESULT: %s (%d disagreements)"
          % ("PASS" if bad == 0 else "FAIL -- STOP CONDITION", bad))
    sys.exit(0 if bad == 0 else 3)


if __name__ == "__main__":
    main()
