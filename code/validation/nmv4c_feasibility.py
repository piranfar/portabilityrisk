"""NM-V4C step 1 -- exposure-only feasibility audit.

Reports SUPPORT for each of the 164 both-context ARG families: how many species, genomes,
BioProjects, occurrences and context blocks carry it, how it is represented in the three
comparison hosts, and which of classes A-E are OBSERVABLE for it.

Deliberately absent, because they are outcome quantities and the design is not frozen yet:
species-specific class proportions, odds ratios, entropy, dominance, direction of effect and
any host-vehicle contrast. This file computes none of them, and a reader can check that from
the code.

Class observability is a boolean support question -- can this family be evaluated in that
architecture at all -- not a proportion.
"""
import argparse, collections, csv, hashlib, json, os, sys

VERSION = "nmv4c_feasibility_v1.0.0"
KLEB = "Klebsiella"
AB = "Acinetobacter baumannii"
PA = "Pseudomonas aeruginosa"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")

    both = [r["gene_family"] for r in csv.DictReader(
        open(os.path.join(O, "both_context_determinants.tsv"), encoding="utf-8"),
        delimiter="\t")]
    print("%s\n  both-context families declared: %d" % (VERSION, len(both)))
    bothset = set(both)

    cls = list(csv.DictReader(open(os.path.join(O, "determinant_portability_classes.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    if len(cls) != 74349:
        print("REFUSING: PRIMARY is %d, expected 74349" % len(cls)); sys.exit(1)
    print("  PRIMARY denominator verified: %d" % len(cls))

    # chromosomal occurrence -> context block, by coordinate containment on its replicon
    blocks_by_rep = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(O, "shared_context_blocks.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        blocks_by_rep[r["replicon_accession"]].append(
            (int(r["block_start"]), int(r["block_end"]), r["block_id"]))
    for v in blocks_by_rep.values():
        v.sort()

    def block_of(rep, s, e):
        for bs, be, bid in blocks_by_rep.get(rep, ()):
            if bs <= s and be >= e:
                return bid
        return None

    F = collections.defaultdict(lambda: {
        "species": set(), "genomes": set(), "bioprojects": set(), "occ": 0,
        "blocks": set(), "classes": set(),
        "ab_occ": 0, "ab_bp": set(), "pa_occ": 0, "pa_bp": set(),
        "kleb_occ": 0, "kleb_bp": set(), "kleb_species": set()})
    unmapped = 0
    for r in cls:
        fam = r["gene_family"]
        if fam not in bothset:
            continue
        d = F[fam]
        sp = r["organism_harmonized"]
        d["species"].add(sp); d["genomes"].add(r["assembly_version"])
        d["bioprojects"].add(r["bioproject_accession"]); d["occ"] += 1
        d["classes"].add(r["portability_class"])
        if r["replicon_molecule_type"].lower().startswith("chrom"):
            b = block_of(r["replicon_accession"], int(r["gene_start"]), int(r["gene_end"]))
            if b:
                d["blocks"].add(b)
            else:
                unmapped += 1
        if sp == AB:
            d["ab_occ"] += 1; d["ab_bp"].add(r["bioproject_accession"])
        elif sp == PA:
            d["pa_occ"] += 1; d["pa_bp"].add(r["bioproject_accession"])
        if r["genus"] == KLEB or sp.startswith(KLEB):
            d["kleb_occ"] += 1; d["kleb_bp"].add(r["bioproject_accession"])
            d["kleb_species"].add(sp)
    print("  chromosomal occurrences with no containing block: %d" % unmapped)
    print("  families resolved from the class table: %d of %d" % (len(F), len(bothset)))

    cols = ["gene_family", "n_species", "n_genomes", "n_bioprojects", "n_occurrences",
            "n_context_blocks", "ab_occurrences", "ab_bioprojects", "pa_occurrences",
            "pa_bioprojects", "kleb_occurrences", "kleb_bioprojects", "kleb_species",
            "class_A_observable", "class_B_observable", "class_C_observable",
            "class_D_observable", "class_E_observable", "n_classes_observable",
            "present_in_ab", "present_in_pa", "present_in_kleb", "present_in_all_three"]
    rows = []
    for fam in sorted(F):
        d = F[fam]
        obs = {c: ("yes" if c in d["classes"] else "no") for c in "ABCDE"}
        rows.append({
            "gene_family": fam, "n_species": len(d["species"]),
            "n_genomes": len(d["genomes"]), "n_bioprojects": len(d["bioprojects"]),
            "n_occurrences": d["occ"], "n_context_blocks": len(d["blocks"]),
            "ab_occurrences": d["ab_occ"], "ab_bioprojects": len(d["ab_bp"]),
            "pa_occurrences": d["pa_occ"], "pa_bioprojects": len(d["pa_bp"]),
            "kleb_occurrences": d["kleb_occ"], "kleb_bioprojects": len(d["kleb_bp"]),
            "kleb_species": len(d["kleb_species"]),
            "class_A_observable": obs["A"], "class_B_observable": obs["B"],
            "class_C_observable": obs["C"], "class_D_observable": obs["D"],
            "class_E_observable": obs["E"], "n_classes_observable": len(d["classes"]),
            "present_in_ab": "yes" if d["ab_occ"] else "no",
            "present_in_pa": "yes" if d["pa_occ"] else "no",
            "present_in_kleb": "yes" if d["kleb_occ"] else "no",
            "present_in_all_three": "yes" if (d["ab_occ"] and d["pa_occ"] and d["kleb_occ"]) else "no"})
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- support landscape, to propose thresholds. No outcome quantity is computed. ----
    print("\n=== SUPPORT LANDSCAPE (exposure only) ===")
    def tally(pred, label):
        n = sum(1 for r in rows if pred(r))
        print("  %-58s %3d of %d" % (label, n, len(rows)))
        return n
    tally(lambda r: r["present_in_all_three"] == "yes", "present in A.b AND P.a AND Klebsiella")
    tally(lambda r: r["present_in_ab"] == "yes", "present in A. baumannii")
    tally(lambda r: r["present_in_pa"] == "yes", "present in P. aeruginosa")
    tally(lambda r: r["present_in_kleb"] == "yes", "present in Klebsiella")
    print()
    for k in (2, 3, 4, 5):
        tally(lambda r, k=k: r["n_species"] >= k, "observed in >=%d species" % k)
    print()
    for k in (5, 10, 20, 30, 50):
        tally(lambda r, k=k: r["n_bioprojects"] >= k, "observed in >=%d BioProjects" % k)
    print()
    for k in (20, 50, 100, 200):
        tally(lambda r, k=k: r["n_occurrences"] >= k, "with >=%d occurrences" % k)
    print()
    tally(lambda r: r["n_classes_observable"] >= 2, "at least 2 of classes A-E observable")
    tally(lambda r: r["n_classes_observable"] >= 3, "at least 3 of classes A-E observable")
    tally(lambda r: r["class_B_observable"] == "yes", "class B observable")
    tally(lambda r: r["class_B_observable"] == "yes"
          and (r["class_C_observable"] == "yes" or r["class_D_observable"] == "yes"
               or r["class_E_observable"] == "yes"),
          "class B AND at least one plasmid class observable")
    print()
    print("  candidate combined support rules:")
    for lab, pred in [
        ("A: >=2 species, >=10 BioProjects, >=20 occurrences",
         lambda r: r["n_species"] >= 2 and r["n_bioprojects"] >= 10 and r["n_occurrences"] >= 20),
        ("B: >=3 species, >=10 BioProjects, >=20 occurrences",
         lambda r: r["n_species"] >= 3 and r["n_bioprojects"] >= 10 and r["n_occurrences"] >= 20),
        ("C: rule B AND class B observable AND >=1 plasmid class observable",
         lambda r: r["n_species"] >= 3 and r["n_bioprojects"] >= 10 and r["n_occurrences"] >= 20
         and r["class_B_observable"] == "yes"
         and any(r["class_%s_observable" % c] == "yes" for c in "CDE")),
        ("D: rule C AND present in A. baumannii",
         lambda r: r["n_species"] >= 3 and r["n_bioprojects"] >= 10 and r["n_occurrences"] >= 20
         and r["class_B_observable"] == "yes"
         and any(r["class_%s_observable" % c] == "yes" for c in "CDE")
         and r["present_in_ab"] == "yes"),
    ]:
        tally(pred, lab)
    print("\n  wrote %s\n  sha256 %s" % (a.out, sha256_file(a.out)))
    print("  NOTE: no class proportion, odds ratio, entropy or contrast was computed here.")


if __name__ == "__main__":
    main()
