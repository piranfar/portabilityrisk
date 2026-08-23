"""NM-V4D freeze -- three-architecture design, fixed before any outcome is scored.

NM-V4C's focused contrast collapsed class A away, and class A is exactly where the two
low-plasmid hosts diverge. That is why gate G5 was qualified (amendment 001) and why this
module exists. The model is widened, not the result reinterpreted after the fact: the three
architectures and all seven success criteria are fixed here, before scoring.
"""
import argparse, csv, datetime, hashlib, json, os, sys

VERSION = "nmv4d_freeze_v1.0.0"
NMV4C_DESIGN_SHA = "8a3c76b157cbf2cd5279342cb2752e1974a2baac976b62455ea0c66f3de4495d"
AMENDMENT = "NMV4C_AMENDMENT_001_COPRIMARY_GATE_CLARIFICATION.md"
AB, PA, KL = "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Klebsiella"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--nmv4c", required=True)
    ap.add_argument("--amendment", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if sha256_file(a.nmv4c) != NMV4C_DESIGN_SHA:
        print("REFUSING: NM-V4C design digest mismatch"); sys.exit(1)
    if os.path.exists(a.out):
        print("REFUSING: %s exists; a frozen design is never overwritten" % a.out); sys.exit(1)
    C = json.load(open(a.nmv4c, encoding="utf-8"))
    ELIG = C["eligibility"]["eligible_families"]
    print("%s | NM-V4C design verified %s" % (VERSION, NMV4C_DESIGN_SHA[:16]))
    print("  inheriting the frozen eligible family set: %d families" % len(ELIG))

    D = {
     "design": "NM_V4D_FROZEN_THREE_ARCHITECTURE_DESIGN",
     "version": "1.0.0", "builder": VERSION,
     "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
     "authorised_by": "owner, 2026-08-21, after NM-V4C acceptance and amendment 001",
     "frozen_before_any_outcome_was_scored": True,
     "supersedes_nothing": "NM-V4C stands unchanged. This module widens the model; it does "
                           "not reinterpret NM-V4C after the fact.",
     "provenance": {"nmv4c_design_sha256": NMV4C_DESIGN_SHA,
                    "amendment": AMENDMENT,
                    "amendment_sha256": sha256_file(a.amendment),
                    "why_this_module_exists":
                      "the NM-V4C focused contrast conditioned away class A, which is the axis "
                      "on which Acinetobacter baumannii and Pseudomonas aeruginosa diverge. "
                      "The two-route model is underspecified for these data."},

     "architectures": {
       "definition": "three mutually exclusive architectures partitioning all 74,349 PRIMARY "
                     "occurrences with nothing discarded",
       "chromosomal_quiescent": {"classes": ["A"],
         "meaning": "chromosomal, no MGE marker within 10 kb on current evidence",
         "caution": "absence of a detected marker is not absence of mobility"},
       "chromosomal_mobile": {"classes": ["B"],
         "meaning": "chromosomal, at least one MGE marker within 10 kb",
         "caution": "proximity is co-location; no transposition or transfer was observed"},
       "plasmid_borne": {"classes": ["C", "D", "E"],
         "meaning": "on a documented plasmid, any mobility-marker status",
         "caution": "class C is retained inside this group and is never dropped; a "
                    "marker-negative plasmid is not a non-mobilizable plasmid"},
       "class_C_retained": True, "A_and_B_never_merged": True,
       "secondary_full_five_class": "the complete A-E composition is reported per family x "
                                    "host and compared by total variation distance"},

     "primary_comparisons": [
       {"id": "PC1", "hosts": [AB, PA], "contrast": ["chromosomal_mobile", "chromosomal_quiescent"],
        "prespecified_direction": "Acinetobacter baumannii enriched for chromosomal_mobile",
        "why": "the comparison NM-V4C could not represent; both hosts are plasmid-poor so the "
               "informative axis is B against A"},
       {"id": "PC2", "hosts": [AB, KL], "contrast": ["chromosomal_mobile", "plasmid_borne"],
        "prespecified_direction": "Acinetobacter baumannii enriched for chromosomal_mobile"},
       {"id": "PC3", "hosts": [PA, KL], "contrast": ["chromosomal_quiescent", "plasmid_borne"],
        "prespecified_direction": "Pseudomonas aeruginosa enriched for chromosomal_quiescent"}],

     "eligibility": {
       "inherited_from_nmv4c": True,
       "rule": ">=3 species AND >=10 BioProjects AND >=20 occurrences AND class B observable "
               "AND at least one of C, D, E observable",
       "per_comparison_additional_rule":
         "a family enters a given comparison only if BOTH hosts have at least one occurrence "
         "in the union of the two architectures being contrasted; otherwise the stratum has a "
         "zero margin and contributes no information",
       "n_eligible": len(ELIG), "eligible_families": ELIG},

     "primary_analysis": {
       "method": "BioProject-balanced, family-stratified Mantel-Haenszel odds ratio with "
                 "Robins-Breslow-Greenland variance",
       "bioproject_balancing_rule":
         "within each (family, host, BioProject) cell the architecture counts are normalised "
         "to unit mass, so every BioProject contributes exactly 1 unit to its (family, host) "
         "stratum regardless of how many genomes it sequenced. Effective sample size is "
         "therefore the NUMBER OF BIOPROJECTS, not the number of occurrences.",
       "why": "an occurrence-level analysis lets one heavily sequenced project inside one host "
              "dominate a family stratum. Balancing deliberately produces WIDER intervals than "
              "the occurrence-level analysis; that is the conservative choice and it is chosen "
              "here before any result is seen.",
       "model_prespecified": "Mantel-Haenszel, the method already established in this project "
                             "for stratified adjustment. No model will be selected on the "
                             "basis of effect strength, and no alternative model will be "
                             "substituted if MH gives a weaker answer."},

     "sensitivity_analyses": [
       {"id": "S1", "name": "occurrence-level", "rule": "raw occurrence counts, unbalanced"},
       {"id": "S2", "name": "one genome per BioProject",
        "rule": "deterministic selection of the assembly with the smallest SHA-256 of its "
                "accession within each BioProject"},
       {"id": "S3", "name": "family-weighted summary",
        "rule": "unweighted mean of per-family log odds ratios, so each family counts once "
                "regardless of its occurrence count; reported beside the MH pooled estimate"},
       {"id": "S4", "name": "leave-one-family-out", "rule": "every eligible family removed once"},
       {"id": "S5", "name": "leave-one-BioProject-out", "rule": "every BioProject removed once"},
       {"id": "S6", "name": "permutation",
        "rule": "host labels permuted over BioProjects within family strata, 2000 permutations, "
                "seed 20260821, two-sided on |ln OR|"}],

     "success_criteria": {
       "SC1": "PC1 odds ratio greater than 1 with a 95 per cent CI excluding 1 - "
              "Acinetobacter baumannii enriched for chromosomal_mobile against "
              "chromosomal_quiescent relative to Pseudomonas aeruginosa",
       "SC2": "PC2 odds ratio greater than 1 with a 95 per cent CI excluding 1 - "
              "Acinetobacter baumannii enriched for chromosomal_mobile against plasmid_borne "
              "relative to Klebsiella",
       "SC3": "PC3 odds ratio greater than 1 with a 95 per cent CI excluding 1 - "
              "Pseudomonas aeruginosa enriched for chromosomal_quiescent against plasmid_borne "
              "relative to Klebsiella",
       "SC4": "all three directions hold in the BioProject-balanced primary AND in the "
              "occurrence-level sensitivity S1 AND in the one-genome-per-BioProject "
              "sensitivity S2",
       "SC5": "no single family or BioProject controls any principal comparison: leave-one-out "
              "relative change in pooled ln odds ratio no greater than 0.20, and the largest "
              "single family carries no more than 30 per cent of the MH weight",
       "SC6": "the full five-class composition is coherent with the three-class result in "
              "direction for all three comparisons",
       "SC7": "an independent verifier reproduces every exported cell and every headline "
              "effect from the exported TSVs",
       "effect_sizes_need_not_match_across_families": True,
       "influence_threshold_relative_log_or": 0.20,
       "max_single_family_weight_share": 0.30},

     "permitted_headline_if_all_criteria_pass":
       "Across matched resistance-gene families, bacterial hosts segregate the same "
       "determinants into three distinct portability architectures: plasmid-conjugative, "
       "chromosomal-mobile and chromosomal-quiescent.",
     "terminology_binding": {
       "preferred": "Identical resistance determinants are routed into distinct "
                    "plasmid-conjugative, chromosomal-mobile and chromosomal-quiescent "
                    "architectures across bacterial hosts.",
       "permitted": "Portability architecture is a host-associated property of the "
                    "determinant-host combination.",
       "host_determines": "usable descriptively ONLY if immediately defined as statistical "
                          "association and not causal control"},
     "prohibited_interpretations": ["observed horizontal transfer", "causal host control",
                                    "demonstrated conjugation", "transfer rate", "clinical risk"],
     "standing_conditionality":
       "every chromosomal_mobile (class B) conclusion remains conditional on NM-V1 validation "
       "of the MGE layer, which has not been run",

     "input_digests": {},
    }
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    for f in ("determinant_portability_classes.tsv", "genome_level_summary.tsv",
              "both_context_determinants.tsv"):
        D["input_digests"][f] = sha256_file(os.path.join(O, f))
    json.dump(D, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("  wrote %s" % a.out)
    print("  FROZEN DESIGN SHA-256: %s" % sha256_file(a.out))


if __name__ == "__main__":
    main()
