"""PR-CONTEXT Phase C -- resistance / stress / biocide / virulence convergence.

Builds FIVE separately labelled evidence layers from the same frozen AMRFinderPlus output and
keeps them apart at every step. A metal-resistance row is never counted as an acquired ARG, a
virulence row is never counted as an acquired ARG, and point mutations are held in their own
layer with their own reduced denominator -- G_C_AMENDMENT_001 records that they were only
searched on the 4,674 genomes that received an organism flag, so their denominator is 2,542
genomes short and pooling them with acquired genes would be a category error.

Co-location on one replicon is reported as CO-LOCATION and nothing more. It is not
co-transfer, it is not co-selection, and it is not a phenotype. No transfer experiment and no
selection experiment exists anywhere in this dataset.
"""
import argparse, collections, csv, glob, hashlib, json, os, sys

VERSION = "pr_context_convergence_v1.0.0"
NEAR_BP = 10000
ADJACENT_BP = 1000


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


def layer_of(t, s, cls, scope=""):
    """The frozen evidence layers. Mutually exclusive by construction.

    CORRECTION applied before any convergence number was reported. The first version mapped
    every non-EFFLUX AMR/AMR row to one 'acquired_amr' layer, which is sensitivity set S1
    (76,383 rows), not the frozen PRIMARY definition (core, non-EFFLUX, 74,349 rows). That
    made the convergence denominator 6,631 ARG-bearing plasmids where every other table in
    this study says 6,621. Two denominators for one quantity is exactly the defect this
    programme keeps having to correct, so the layers now separate scope explicitly and the
    convergence analysis runs on PRIMARY, with S1 reported beside it as a sensitivity.
    """
    if t == "AMR" and s == "AMR":
        if cls == "EFFLUX":
            return "acquired_amr_efflux"
        return "acquired_amr" if scope == "core" else "acquired_amr_plus_nonefflux"
    if t == "AMR" and s in ("POINT", "POINT_DISRUPT"):
        return "point_mutation"
    if t == "STRESS" and s == "METAL":
        return "metal_stress"
    if t == "STRESS" and s == "BIOCIDE":
        return "biocide"
    if t == "STRESS":
        return "other_stress"
    if t == "VIRULENCE":
        return "virulence"
    return "other"


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
    mob = {r["replicon_accession"]: r for r in csv.DictReader(
        open(os.path.join(out, "plasmid_mobility_annotation.tsv"), encoding="utf-8"),
        delimiter="\t")}

    feats = collections.defaultdict(list)      # (asm, replicon) -> feature list
    layer_rows = collections.Counter()
    layer_by_mol = collections.Counter()
    for f in sorted(glob.glob(os.path.join(a.amrfinder_dir, "*.tsv"))):
        if f.endswith("_timing.tsv"):
            continue
        asm = os.path.basename(f)[:-4]
        idx = inv.get(asm, {})
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"),
                                delimiter="\t"):
            L = layer_of(r.get("Type"), r.get("Subtype"), r.get("Class"),
                         r.get("Scope", ""))
            layer_rows[L] += 1
            cid = (r.get("Contig id") or "").strip()
            rec = idx.get(cid) or idx.get(cid.split(".")[0])
            mt = (rec or {}).get("replicon_molecule_type", "") or "unresolved"
            layer_by_mol[(L, mt)] += 1
            if rec is None:
                continue
            try:
                s, e = sorted((int(r["Start"]), int(r["Stop"])))
            except Exception:
                continue
            feats[(asm, rec["sequence_accession"])].append({
                "layer": L, "symbol": r.get("Element symbol", ""),
                "cls": r.get("Class", ""), "subclass": r.get("Subclass", ""),
                "start": s, "end": e, "mol": mt,
                "genus": "", "organism": ""})

    # attach organism from the occurrence table
    org = {}
    for r in csv.DictReader(open(os.path.join(out, "determinant_occurrences.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        org[r["assembly_version"]] = (r["genus"], r["organism_harmonized"],
                                      r["biosample_accession"],
                                      r["bioproject_accession"])

    lrows = [{"evidence_layer": k, "n_rows": v,
              "counted_as_acquired_arg": "yes" if k == "acquired_amr" else "NO",
              "note": {"acquired_amr": "the frozen PRIMARY layer (core, non-EFFLUX) - "
                                       "the convergence denominator",
                       "acquired_amr_plus_nonefflux": "plus-scope non-EFFLUX; sensitivity "
                                                      "set S1 only, NOT in the primary",
                       "acquired_amr_efflux": "EFFLUX class; sensitivity set S2 only",
                       "point_mutation": "reduced denominator - only searched where an "
                                         "AMRFinderPlus organism flag existed",
                       "metal_stress": "metal tolerance, NOT antibiotic resistance",
                       "biocide": "biocide tolerance, NOT antibiotic resistance",
                       "other_stress": "heat/acid tolerance, NOT antibiotic resistance",
                       "virulence": "virulence, NOT antibiotic resistance"}.get(k, "")}
             for k, v in sorted(layer_rows.items(), key=lambda kv: -kv[1])]
    write(os.path.join(out, "evidence_layer_inventory.tsv"), lrows, list(lrows[0].keys()))

    # ---------------- replicon-level convergence
    prow = []
    for (asm, acc), fl in sorted(feats.items()):
        mol = fl[0]["mol"]
        if mol != "Plasmid":
            continue
        arg = [x for x in fl if x["layer"] == "acquired_amr"]
        if not arg:
            continue
        met = [x for x in fl if x["layer"] == "metal_stress"]
        bio = [x for x in fl if x["layer"] == "biocide"]
        vir = [x for x in fl if x["layer"] == "virulence"]
        classes = {x["cls"] for x in arg if x["cls"]}
        m = mob.get(acc, {})
        cat = m.get("portability_mobility_category", "mobility_unresolved")

        def prox(a_list, b_list):
            """closest relationship between two feature sets on the same replicon"""
            if not a_list or not b_list:
                return ""
            best = None
            for x in a_list:
                for y in b_list:
                    gap = max(0, max(x["start"], y["start"]) - min(x["end"], y["end"]))
                    if best is None or gap < best:
                        best = gap
            if best is None:
                return ""
            if best == 0:
                return "overlapping_or_abutting"
            if best <= ADJACENT_BP:
                return "adjacent_within_1kb"
            if best <= NEAR_BP:
                return "nearby_within_10kb"
            return "same_replicon_only"
        g, o, bs, bp = org.get(asm, ("", "", "", ""))
        prow.append({
            "assembly_version": asm, "replicon_accession": acc,
            "biosample_accession": bs, "bioproject_accession": bp,
            "organism_harmonized": o, "genus": g,
            "plasmid_length": m.get("plasmid_length", ""),
            "rep_types": m.get("rep_types", ""),
            "mobility_category": cat,
            "n_acquired_arg": len(arg),
            "n_arg_drug_classes": len(classes),
            "arg_drug_classes": ";".join(sorted(classes)),
            "arg_symbols": ";".join(sorted({x["symbol"] for x in arg})),
            "n_metal_stress": len(met), "n_biocide": len(bio), "n_virulence": len(vir),
            "multi_arg_ge2": "yes" if len(arg) >= 2 else "no",
            "multi_class_ge2": "yes" if len(classes) >= 2 else "no",
            "multi_class_ge3": "yes" if len(classes) >= 3 else "no",
            "arg_plus_metal": "yes" if met else "no",
            "arg_plus_biocide": "yes" if bio else "no",
            "arg_plus_virulence": "yes" if vir else "no",
            "arg_stress_virulence_triple":
                "yes" if (met or bio) and vir else "no",
            "proximity_arg_metal": prox(arg, met),
            "proximity_arg_biocide": prox(arg, bio),
            "proximity_arg_virulence": prox(arg, vir),
            "evidence_statement":
                "CO-LOCATION on one documented replicon. Not co-transfer, not co-selection, "
                "not a phenotype."})
    write(os.path.join(out, "plasmid_convergence.tsv"), prow, list(prow[0].keys()))

    n = len(prow)
    def pct(k, v="yes"):
        c = sum(1 for r in prow if r[k] == v)
        return c, 100.0 * c / n
    summ = []
    for lab, k in (("multiple acquired ARGs (>=2)", "multi_arg_ge2"),
                   ("ARGs from >=2 drug classes", "multi_class_ge2"),
                   ("ARGs from >=3 drug classes", "multi_class_ge3"),
                   ("ARG + metal-resistance marker", "arg_plus_metal"),
                   ("ARG + biocide marker", "arg_plus_biocide"),
                   ("ARG + virulence marker", "arg_plus_virulence"),
                   ("ARG + stress + virulence", "arg_stress_virulence_triple")):
        c, p = pct(k)
        summ.append({"convergence_feature": lab, "n_plasmids": c,
                     "denominator_unique_arg_plasmids": n, "pct": round(p, 2)})
    s1 = set()
    for (asm, acc), fl in feats.items():
        if fl[0]["mol"] == "Plasmid" and any(
                x["layer"] in ("acquired_amr", "acquired_amr_plus_nonefflux") for x in fl):
            s1.add(acc)
    summ.append({"convergence_feature":
                 "SENSITIVITY S1: unique ARG-bearing plasmids if plus-scope non-EFFLUX "
                 "genes are also admitted",
                 "n_plasmids": len(s1), "denominator_unique_arg_plasmids": n,
                 "pct": ""})
    write(os.path.join(out, "plasmid_convergence_summary.tsv"), summ, list(summ[0].keys()))

    # by mobility category
    byc = []
    for cat in ("predicted_conjugative", "predicted_mobilizable",
                "nonconjugative_or_no_mobility_markers_detected"):
        sub = [r for r in prow if r["mobility_category"] == cat]
        if not sub:
            continue
        d = len(sub)
        row = {"mobility_category": cat, "n_plasmids": d}
        for lab, k in (("multi_arg_ge2", "multi_arg_ge2"),
                       ("multi_class_ge2", "multi_class_ge2"),
                       ("multi_class_ge3", "multi_class_ge3"),
                       ("arg_plus_metal", "arg_plus_metal"),
                       ("arg_plus_biocide", "arg_plus_biocide"),
                       ("arg_plus_virulence", "arg_plus_virulence"),
                       ("arg_stress_virulence_triple", "arg_stress_virulence_triple")):
            c = sum(1 for r in sub if r[k] == "yes")
            row["pct_" + lab] = round(100.0 * c / d, 2)
            row["n_" + lab] = c
        row["median_arg_per_plasmid"] = sorted(
            int(r["n_acquired_arg"]) for r in sub)[d // 2]
        byc.append(row)
    write(os.path.join(out, "convergence_by_mobility_class.tsv"), byc,
          list(byc[0].keys()))

    # high-concern architectures
    hc = [r for r in prow
          if r["mobility_category"] == "predicted_conjugative"
          and r["multi_class_ge3"] == "yes"
          and (r["arg_plus_metal"] == "yes" or r["arg_plus_virulence"] == "yes")]
    arch = collections.Counter()
    for r in hc:
        key = (r["rep_types"].split(",")[0].strip() if r["rep_types"] else "untyped",
               r["arg_drug_classes"])
        arch[key] += 1
    arows = [{"rep_type": k[0], "drug_class_combination": k[1], "n_plasmids": v,
              "criteria": "predicted conjugative AND >=3 drug classes AND (metal or "
                          "virulence co-location)",
              "evidence_type": "documented replicon location + database-predicted mobility "
                               "+ co-located sequence features",
              "not_claimed": "co-transfer, co-selection, or demonstrated phenotype"}
             for k, v in arch.most_common(40)]
    write(os.path.join(out, "high_concern_plasmid_architectures.tsv"), arows,
          list(arows[0].keys()) if arows else ["rep_type"])

    print("%s\n=== PHASE C: five separated evidence layers ===" % VERSION)
    for r in lrows:
        print("  %-22s %7d rows   counted as acquired ARG: %s"
              % (r["evidence_layer"], r["n_rows"], r["counted_as_acquired_arg"]))
    print("\n=== convergence over %d unique ARG-bearing plasmids ===" % n)
    for r in summ:
        print("  %-34s %5d  %s" % (r["convergence_feature"][:34], r["n_plasmids"],
                                   ("%6.2f%%" % r["pct"]) if r["pct"] != "" else "(sensitivity)"))
    print("\n=== convergence by predicted mobility ===")
    for r in byc:
        print("  %-48s n=%-5d  >=3 classes %5.1f%%  +metal %5.1f%%  +vir %5.1f%%  "
              "triple %5.1f%%  median ARG %d"
              % (r["mobility_category"], r["n_plasmids"], r["pct_multi_class_ge3"],
                 r["pct_arg_plus_metal"], r["pct_arg_plus_virulence"],
                 r["pct_arg_stress_virulence_triple"], r["median_arg_per_plasmid"]))
    print("\nhigh-concern architectures: %d plasmids, %d distinct rep-type/class "
          "combinations" % (len(hc), len(arch)))
    for r in arows[:8]:
        print("  %-16s %-58s %d" % (r["rep_type"], r["drug_class_combination"][:58],
                                    r["n_plasmids"]))
    for f in ("evidence_layer_inventory.tsv", "plasmid_convergence.tsv",
              "plasmid_convergence_summary.tsv", "convergence_by_mobility_class.tsv",
              "high_concern_plasmid_architectures.tsv"):
        print("  %-46s %s" % (f, sha256_file(os.path.join(out, f))))


if __name__ == "__main__":
    main()
