"""PR-CONTEXT step 3 -- the direct determinant-to-replicon join.

Every acquired ARG occurrence is matched to a DOCUMENTED replicon record by exact accession.
Nothing is predicted, nothing is inferred from a gene name, and no occurrence is dropped: a
row that cannot be matched keeps its place in the missingness denominator carrying the reason
it could not be matched.

The rule that matters most here is the one that is easiest to get wrong: absence of a plasmid
match is NOT evidence of chromosomal location. An occurrence whose contig id resolves to no
replicon record is identifier_unmatched, not chromosomal.
"""
import argparse, collections, csv, glob, hashlib, json, os, re, sys, datetime

VERSION = "pr_context_build_join_v1.0.0"

# Frozen in AMENDMENT_001, written before any enrichment statistic was computed.
_ALLELE_TAIL = re.compile(r"^(?:\d+[A-Za-z'\"]*|[IVX]+[a-z]?|[A-Za-z]{1,3})$")


def gene_family(symbol):
    """Strip a terminal ALLELE designator from an AMRFinderPlus element symbol.

    AMRFinderPlus 4.2.7 nucleotide output carries no family column, and the run may not be
    repeated to add one, so the family is derived by a fixed textual rule and the raw symbol
    is always retained beside it.

    blaTEM-1 -> blaTEM ; blaCTX-M-15 -> blaCTX-M ; aph(6)-Id -> aph(6) ;
    aph(3'')-Ib -> aph(3'') ; sul1 -> sul1 ; tet(A) -> tet(A)
    """
    s = (symbol or "").strip()
    if "-" not in s:
        return s
    head, _, tail = s.rpartition("-")
    if head and _ALLELE_TAIL.match(tail):
        return head
    return s


AMR_COLS = ["assembly_accession", "assembly_version", "biosample_accession",
            "bioproject_accession", "organism_original", "organism_harmonized", "genus",
            "sequence_accession", "replicon_accession", "replicon_name", "replicon_length",
            "replicon_molecule_type", "determinant_name", "gene_family", "drug_class",
            "subclass", "scope", "gene_start", "gene_end", "strand", "identity", "coverage",
            "amrfinder_method", "element_name", "call_completeness", "analysis_set",
            "source_file_sha256", "join_method", "evidence_type", "evidence_status"]


def write(path, rows, cols):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in cols) + "\n")


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--amrfinder-dir", required=True)
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--protocol-sha256", required=True)
    a = ap.parse_args()
    out = os.path.join(a.root, "out")
    pp = os.path.join(out, "FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1.json")
    if sha256_file(pp) != a.protocol_sha256:
        print("REFUSING: frozen protocol digest mismatch"); sys.exit(1)
    P = json.load(open(pp, encoding="utf-8"))
    print("%s | protocol %s verified" % (VERSION, a.protocol_sha256[:16]))

    wl = {r["queried_accession"]: r for r in csv.DictReader(
        open(a.worklist, encoding="utf-8"), delimiter="\t")}
    # replicon index: (assembly, accession-without-version) -> list of records
    rep = collections.defaultdict(list)
    rep_rows = 0
    for r in csv.DictReader(open(os.path.join(out, "replicon_inventory.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        rep_rows += 1
        for key in (r.get("sequence_accession"), r.get("refseq_accession"),
                    r.get("genbank_accession")):
            if key:
                rep[(r["assembly_accession"], key)].append(r)
                rep[(r["assembly_accession"], key.split(".")[0])].append(r)
    print("replicon inventory rows: %d | index keys: %d" % (rep_rows, len(rep)))

    sidecars = {}
    for f in glob.glob(os.path.join(a.amrfinder_dir, "*.run.json")):
        d = json.load(open(f, encoding="utf-8"))
        sidecars[d["queried_accession"]] = d

    files = [f for f in glob.glob(os.path.join(a.amrfinder_dir, "*.tsv"))
             if not f.endswith("_timing.tsv")]
    rows, audit = [], collections.Counter()
    n_all_amr = 0
    for f in sorted(files):
        acc = os.path.basename(f)[:-4]
        w = wl.get(acc, {})
        sc = sidecars.get(acc, {})
        sp = (w.get("species") or sc.get("species") or "").strip()
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"), delimiter="\t"):
            if r.get("Type") != "AMR" or r.get("Subtype") != "AMR":
                continue
            n_all_amr += 1
            cls = r.get("Class") or ""
            scope = r.get("Scope") or ""
            if scope == "core" and cls != "EFFLUX":
                aset = "PRIMARY"
            elif cls != "EFFLUX":
                aset = "S1_plus_nonefflux"
            else:
                aset = "S2_efflux"
            cid = (r.get("Contig id") or "").strip()
            cands = rep.get((acc, cid)) or rep.get((acc, cid.split(".")[0])) or []
            uniq = {c["sequence_accession"]: c for c in cands}
            method = r.get("Method", "")
            try:
                gs, ge = int(r["Start"]), int(r["Stop"])
                coords_ok = True
            except Exception:
                gs = ge = ""
                coords_ok = False

            if not coords_ok:
                ev, st, jm, rec = "coordinates_missing", "UNRESOLVED", "none", {}
            elif not uniq:
                ev, st, jm, rec = "identifier_unmatched", "UNRESOLVED", "none", {}
            elif len(uniq) > 1:
                ev, st, jm, rec = "identifier_ambiguous", "UNRESOLVED", "multi", {}
            else:
                rec = list(uniq.values())[0]
                jm = "exact_accession_version" if cid == rec.get("sequence_accession") \
                    else "exact_accession_unversioned"
                mt = (rec.get("replicon_molecule_type") or "").strip()
                L = rec.get("replicon_length")
                try:
                    within = 1 <= min(gs, ge) and max(gs, ge) <= int(L)
                except Exception:
                    within = False
                if not within:
                    ev, st = "coordinates_missing", "UNRESOLVED"
                    jm = jm + "|coords_outside_replicon"
                elif mt == "Chromosome":
                    ev, st = "direct_chromosome", "RESOLVED"
                elif mt == "Plasmid":
                    ev, st = "direct_plasmid", "RESOLVED"
                elif mt:
                    ev, st = "direct_other", "RESOLVED_OTHER"
                else:
                    ev, st = "replicon_unclassified", "UNRESOLVED"
            audit[(aset, ev)] += 1
            rows.append({
                "assembly_accession": acc.split(".")[0], "assembly_version": acc,
                "biosample_accession": w.get("biosample_accession", ""),
                "bioproject_accession": w.get("bioproject_accession", ""),
                "organism_original": sp,
                "organism_harmonized": " ".join(sp.split()[:2]),
                "genus": sp.split()[0] if sp else "",
                "sequence_accession": cid,
                "replicon_accession": rec.get("sequence_accession", ""),
                "replicon_name": rec.get("replicon_name", ""),
                "replicon_length": rec.get("replicon_length", ""),
                "replicon_molecule_type": rec.get("replicon_molecule_type", ""),
                "determinant_name": r.get("Element symbol", ""),
                "gene_family": gene_family(r.get("Element symbol", "")),
                "drug_class": cls, "subclass": r.get("Subclass", ""), "scope": scope,
                "gene_start": gs, "gene_end": ge, "strand": r.get("Strand", ""),
                "identity": r.get("% Identity to reference", ""),
                "coverage": r.get("% Coverage of reference", ""),
                "amrfinder_method": method, "element_name": r.get("Element name", ""),
                "call_completeness": ("PARTIAL" if "PARTIAL" in method else
                                      "INTERNAL_STOP" if method == "INTERNAL_STOP"
                                      else "complete"),
                "analysis_set": aset,
                "source_file_sha256": sc.get("output_sha256", ""),
                "join_method": jm, "evidence_type": ev, "evidence_status": st})

    write(os.path.join(out, "determinant_occurrences.tsv"), rows, AMR_COLS)
    unres = [r for r in rows if r["evidence_status"] != "RESOLVED"]
    write(os.path.join(out, "unmatched_or_ambiguous_determinants.tsv"), unres, AMR_COLS)

    prim = [r for r in rows if r["analysis_set"] == "PRIMARY"]
    print("\nAMR/AMR rows read: %d" % n_all_amr)
    print("PRIMARY set (core, non-EFFLUX): %d  (protocol declared %d)"
          % (len(prim), P["acquired_arg_definition"]["primary_n_expected"]))
    if len(prim) != P["acquired_arg_definition"]["primary_n_expected"]:
        print("STOP CONDITION: primary denominator differs from the frozen declaration")
        sys.exit(2)

    ac = collections.Counter(r["evidence_type"] for r in prim)
    print("\n=== Q1 PRIMARY: direct replicon evidence over %d acquired ARG occurrences ==="
          % len(prim))
    for k in ("direct_chromosome", "direct_plasmid", "direct_other",
              "replicon_unclassified", "identifier_ambiguous", "identifier_unmatched",
              "coordinates_missing"):
        n = ac.get(k, 0)
        print("  %-24s %7d  %6.2f%%" % (k, n, 100.0 * n / len(prim)))
    res = ac.get("direct_chromosome", 0) + ac.get("direct_plasmid", 0)
    if res:
        print("  --- of the %d with unambiguous chromosome/plasmid origin:" % res)
        print("      plasmid    %7d  %6.2f%%" % (ac.get("direct_plasmid", 0),
                                                 100.0 * ac.get("direct_plasmid", 0) / res))
        print("      chromosome %7d  %6.2f%%" % (ac.get("direct_chromosome", 0),
                                                 100.0 * ac.get("direct_chromosome", 0) / res))

    arows = [{"analysis_set": k[0], "evidence_type": k[1], "n": v} for k, v in
             sorted(audit.items())]
    write(os.path.join(out, "determinant_replicon_join_audit.tsv"), arows,
          ["analysis_set", "evidence_type", "n"])
    cov = []
    for aset in ("PRIMARY", "S1_plus_nonefflux", "S2_efflux"):
        sub = [r for r in rows if r["analysis_set"] == aset]
        c = collections.Counter(r["evidence_type"] for r in sub)
        r_ = c.get("direct_chromosome", 0) + c.get("direct_plasmid", 0)
        cov.append({"analysis_set": aset, "n_occurrences": len(sub),
                    "n_resolved_chrom_or_plasmid": r_,
                    "pct_resolved": round(100.0 * r_ / max(len(sub), 1), 4),
                    "n_plasmid": c.get("direct_plasmid", 0),
                    "n_chromosome": c.get("direct_chromosome", 0),
                    "n_unresolved": len(sub) - r_ - c.get("direct_other", 0),
                    "n_other_molecule": c.get("direct_other", 0)})
    write(os.path.join(out, "replicon_mapping_coverage.tsv"), cov, list(cov[0].keys()))

    for f in ("determinant_occurrences.tsv", "determinant_replicon_join_audit.tsv",
              "replicon_mapping_coverage.tsv", "unmatched_or_ambiguous_determinants.tsv"):
        print("  %-46s %s" % (f, sha256_file(os.path.join(out, f))))


if __name__ == "__main__":
    main()
