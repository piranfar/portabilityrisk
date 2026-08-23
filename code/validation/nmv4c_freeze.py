"""NM-V4C freeze -- fix the design before any host-vehicle outcome is scored.

Thresholds are justified by provenance, not by effect size. The occurrence floor of 20 is the
project's OWN pre-existing rule: determinant_plasmid_enrichment.tsv already declines to test
families below 20 occurrences. Reusing a threshold that was frozen for a different purpose,
before this module existed, is the strongest available guarantee that it was not tuned here.

The feasibility audit that informed these thresholds computed only exposure counts and class
OBSERVABILITY booleans. It computed no proportion, odds ratio, entropy or contrast.
"""
import argparse, csv, datetime, hashlib, json, os, sys

VERSION = "nmv4c_freeze_v1.0.0"
PROTOCOL_SHA = "c968fb6d16a528a64d064d6f8bbac745390804df4ad0897c09b12de84ca3fbff"
AUDIT_SHA = "6d6ddf4b6d66cf743470b0625a6ebf266d893cc849d8664d567b7f509754fc13"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--audit", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if sha256_file(a.protocol) != PROTOCOL_SHA:
        print("REFUSING: protocol digest mismatch"); sys.exit(1)
    if sha256_file(a.audit) != AUDIT_SHA:
        print("REFUSING: feasibility audit digest mismatch"); sys.exit(1)
    if os.path.exists(a.out):
        print("REFUSING: %s exists; a frozen design is never overwritten" % a.out); sys.exit(1)
    print("%s | protocol %s | audit %s verified"
          % (VERSION, PROTOCOL_SHA[:12], AUDIT_SHA[:12]))

    rows = list(csv.DictReader(open(a.audit, encoding="utf-8"), delimiter="\t"))

    def elig(r):
        return (int(r["n_species"]) >= 3 and int(r["n_bioprojects"]) >= 10
                and int(r["n_occurrences"]) >= 20
                and r["class_B_observable"] == "yes"
                and any(r["class_%s_observable" % c] == "yes" for c in "CDE"))

    primary = sorted(r["gene_family"] for r in rows if elig(r))
    ab_set = sorted(r["gene_family"] for r in rows if elig(r) and r["present_in_ab"] == "yes")
    print("  eligible families (primary)          : %d" % len(primary))
    print("  eligible AND present in A. baumannii : %d" % len(ab_set))

    D = {
     "design": "NM_V4C_FROZEN_DESIGN",
     "version": "1.0.0",
     "builder": VERSION,
     "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
     "authorised_by": "owner, 2026-08-21, after D-NM10 provisionally resolved for Nature Microbiology",
     "protocol_sha256": PROTOCOL_SHA,
     "feasibility_audit_sha256": AUDIT_SHA,
     "frozen_before_any_outcome_was_scored": True,

     "question": ("The same ARG family can occupy different portability architectures in "
                  "different bacterial hosts; therefore species differences are not explained "
                  "solely by differences in ARG-family composition."),
     "conceptual_source": ("Transferable structure from Mulkern et al., Nat Commun, "
                           "doi:10.1038/s41467-026-76357-y: the determinant alone does not "
                           "specify the biological route; the determinant-host-vehicle "
                           "combination does. No data, result or method is imported from that "
                           "paper; only the analytical framing."),

     "eligibility": {
       "rule": ">=3 species AND >=10 BioProjects AND >=20 occurrences AND class B observable "
               "AND at least one of classes C, D, E observable",
       "justification": {
         "occurrence_floor_20": "the project's own pre-existing rule - families below 20 "
                                "occurrences are not tested in determinant_plasmid_enrichment.tsv. "
                                "Frozen before this module existed and reused unchanged.",
         "species_floor_3": "a host-vehicle contrast is undefined with fewer than 3 hosts; 2 "
                            "hosts permit a difference but no heterogeneity assessment",
         "bioproject_floor_10": "the BioProject cluster bootstrap and the BioProject-balanced "
                                "aggregation both need enough clusters to be meaningful",
         "class_observability": "the focused route contrast B against C+D+E is undefined for a "
                                "family in which neither side is observed"},
       "not_chosen_to_maximise_an_effect": "no proportion, odds ratio or contrast was computed "
                                           "before these thresholds were fixed",
       "n_eligible_primary": len(primary),
       "n_eligible_present_in_ab": len(ab_set),
       "eligible_families": primary,
       "eligible_families_present_in_ab": ab_set},

     "primary_unit": {
       "unit": "ARG-family x species x BioProject",
       "not": "raw occurrence",
       "why": "conditioning on ARG family removes composition differences between hosts, which "
              "is the whole point of the question; conditioning on BioProject removes "
              "pseudo-replication from genomes submitted together. A raw-occurrence unit would "
              "let one heavily sequenced project inside one host dominate a family.",
       "co_primary_operationalisation": [
         "P1-full: Mantel-Haenszel over family strata using all occurrences",
         "P2-balanced: the identical analysis after collapsing to one genome per BioProject "
         "(deterministic: smallest SHA-256 of the assembly accession), so each "
         "family x species x BioProject cell contributes at most one genome"],
       "both_must_agree": "P1 and P2 are CO-PRIMARY. If they disagree in direction or "
                          "significance, the result is reported as sampling-structure "
                          "dependent and the universal claim is not made."},

     "architecture": {
       "five_classes_retained": {"A": "chromosome, no nearby MGE marker",
                                 "B": "chromosome with nearby MGE marker",
                                 "C": "plasmid, marker-negative",
                                 "D": "mobilization-consistent plasmid",
                                 "E": "conjugation-consistent plasmid"},
       "no_class_is_discarded": True,
       "analysis_1_full_five_class": "per family x species, the full A-E composition is "
                                     "reported and compared across hosts by total variation "
                                     "distance; class A is never dropped",
       "analysis_2_focused_route_contrast": {
         "contrast": "B against C+D+E",
         "class_A_handling": "class A is EXCLUDED from the 2x2 but is reported separately as a "
                             "competing context for every family x species cell, with its own "
                             "counts and fraction. It is not silently discarded: a host shift "
                             "from A to B is a different biological statement from a shift "
                             "from plasmid to B, and both are reported.",
         "why_this_contrast": "it is the route question - given that a determinant is in this "
                              "host, does it sit in chromosomal mobile context or on a plasmid"}},

     "primary_test": {
       "method": "family-stratified Mantel-Haenszel odds ratio with Robins-Breslow-Greenland "
                 "variance, strata = ARG family, exposure = host species, outcome = class B "
                 "versus classes C+D+E",
       "why_this_method": "it is the method this project already uses for species-adjusted "
                          "enrichment (determinant_enrichment_species_adjusted.tsv). Reusing "
                          "an established in-project method, with family and species swapped "
                          "in their roles, avoids choosing a method for its result.",
       "host_contrasts": [["Acinetobacter baumannii", "Klebsiella"],
                          ["Acinetobacter baumannii", "Pseudomonas aeruginosa"],
                          ["Pseudomonas aeruginosa", "Klebsiella"]],
       "klebsiella_definition": "genus Klebsiella, all species pooled, because the comparison "
                                "is host-genus level and pooling is declared in advance",
       "global_test": {"method": "permutation of species labels within (family, BioProject) "
                                 "strata, recomputing the pooled MH statistic",
                       "n_permutations": 2000, "seed": 20260821,
                       "why": "defensible for sparse clustered compositional data and makes no "
                              "distributional assumption"},
       "sensitivity_methods": [
         "family fixed-effect logistic model, reported only as a check on the MH direction",
         "full five-class total-variation distance between hosts within family",
         "restriction to families present in all three comparison hosts"],
       "method_was_not_selected_on_results": True},

     "heterogeneity": {
       "required": "quantified, not suppressed",
       "statistics": ["per-family odds ratio with CI", "count of families in each direction",
                      "Cochran Q and I-squared across family strata",
                      "share of total MH weight carried by the single largest family"],
       "note": "families are NOT required to move in the same direction. Heterogeneity is "
               "biologically expected: a family whose route is fixed by its own biology should "
               "not shift with host."},

     "ab_specific_test": {
       "hypothesis": "the same ARG families found in other hosts shift, in A. baumannii, "
                     "towards class B and away from plasmid classes C-E",
       "direction_prespecified": "MH odds ratio for (B vs C+D+E) with A. baumannii as the "
                                 "exposed host is predicted to be GREATER than 1",
       "required_controls": [
         "ARG-family conditioning - the MH strata",
         "BioProject-balanced analysis - co-primary P2",
         "one-genome-per-BioProject sensitivity - identical to P2",
         "exclusion of the most abundant family, one at a time",
         "leave-one-family-out influence over every eligible family",
         "comparison with low-plasmid Pseudomonas aeruginosa",
         "comparison with plasmid-rich Klebsiella"],
       "minimum_supporting_families": 10},

     "gates": {
       "G1_pooled_effect": "the pooled family-conditioned MH odds ratio for A. baumannii "
                           "against Klebsiella has a 95 per cent CI excluding 1",
       "G2_multiple_families": "at least 10 eligible families with an estimable odds ratio "
                               "point in the predicted direction, AND the single largest "
                               "family carries no more than 30 per cent of the MH weight",
       "G3_bioproject_balancing": "the co-primary P2 analysis keeps the same direction with a "
                                  "95 per cent CI excluding 1",
       "G4_leave_one_family_out": "removing any single family moves pooled log odds ratio by "
                                  "no more than 20 per cent relative",
       "G5_pa_control": "the MH odds ratio for A. baumannii against Pseudomonas aeruginosa "
                        "excludes 1 in the same direction, so the effect is not a generic "
                        "low-plasmid artefact",
       "G6_representation_coherence": "the five-class and focused-route representations agree "
                                      "in direction for A. baumannii",
       "G7_bioproject_influence": "removing any single BioProject moves pooled log odds ratio "
                                  "by no more than 20 per cent relative",
       "influence_threshold_relative_log_or": 0.20,
       "all_thresholds_frozen_before_scoring": True},

     "secondary_classifications": {
       "note": "secondary; cannot replace the primary family-conditioned test",
       "computed_on": "families eligible and present in at least 2 comparison hosts",
       "vehicle_stable": "maximum pairwise absolute difference in class-B fraction across "
                         "comparison hosts is below 0.15",
       "host_dependent": "maximum pairwise absolute difference in class-B fraction across "
                         "comparison hosts is at least 0.30, with at least 2 BioProjects "
                         "supporting each host in the pair",
       "broad_route_portable": "at least 4 of the 5 classes observable AND observed in at "
                               "least 5 species",
       "lineage_or_project_private": "at least 70 per cent of the family's occurrences come "
                                     "from a single BioProject",
       "precedence": "lineage_or_project_private is assigned first and overrides the others, "
                     "because a project-private family cannot support a host statement",
       "thresholds_frozen_before_outcome_inspection": True},

     "prohibited_interpretations": [
       "observed horizontal transfer", "causal host control",
       "demonstrated conjugation", "transfer rate", "clinical risk",
       "any claim that a host determines where a gene goes"],
     "permitted_interpretation_if_successful":
       "Portability architecture is a property of the determinant-host combination: identical "
       "ARG families occupy systematically different chromosomal-MGE and plasmid mobility "
       "contexts across bacterial hosts.",
     "if_pooled_fails_but_families_heterogeneous":
       "report host-dependent families as exploratory and do NOT make the universal claim",

     "input_digests": {},
    }
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    for f in ("determinant_portability_classes.tsv", "both_context_determinants.tsv",
              "genome_level_summary.tsv", "shared_context_blocks.tsv",
              "mge_feature_inventory.tsv", "replicon_level_summary.tsv"):
        D["input_digests"][f] = sha256_file(os.path.join(O, f))

    json.dump(D, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  wrote %s" % a.out)
    print("  FROZEN DESIGN SHA-256: %s" % sha256_file(a.out))


if __name__ == "__main__":
    main()
