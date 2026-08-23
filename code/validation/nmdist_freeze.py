"""NM-DIST design freeze -- written before any species-specific distance outcome is computed.

Phase 0 established feasibility and reproduced the frozen nearest_mge_distance_bp column exactly
(16,303 positive / 18,837 censored, 0 mismatches). Only pooled reconciliation counts were seen;
no group-specific distance distribution was calculated.
"""
import argparse, collections, csv, datetime, hashlib, json, os, sys

VERSION = "nmdist_freeze_v1.0.0"
SEED = 20260822

INPUTS = [
 ("audit/data/derived/pr_context/out/arg_neighbourhood_windows.tsv",
  "occurrence to block key; ARG coordinates, replicon length, topology, wrap and truncation flags"),
 ("audit/data/derived/pr_context/out/mge_feature_inventory.tsv",
  "frozen MGE feature calls with coordinates; the endpoint source, NOT redefined here"),
 ("audit/data/derived/pr_context/out/shared_context_blocks.tsv", "the 21,955 frozen blocks"),
 ("audit/data/derived/pr_context/out/arg_mge_neighbourhood.tsv",
  "prior per-occurrence nearest-distance column, used only to validate the convention"),
 ("audit/data/derived/pr_context/out/determinant_portability_classes.tsv",
  "frozen cohort; supplies species, BioProject and class A/B membership"),
]


def sha(p):
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
    if os.path.exists(a.out):
        print("REFUSING: %s exists; a freeze never overwrites one" % a.out); sys.exit(1)

    dig = {}
    for rel, role in INPUTS:
        p = os.path.join(a.repo, rel)
        if not os.path.exists(p):
            print("REFUSING: missing %s" % rel); sys.exit(1)
        dig[rel] = {"sha256": sha(p), "role": role}

    Q = {
 "protocol": "NMDIST_DISTANCE_TO_MGE", "version": "1.0.0", "builder": VERSION,
 "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
 "frozen_before_any_species_specific_distance_outcome": True,

 "purpose": "Replace reliance on a single binary +/-10 kb threshold with a full spatial "
            "description of ARG-to-MGE distance. This is a STRENGTHENING analysis. It does not "
            "alter the frozen PRIMARY cohort, the portability classes or any existing NM result.",

 "PHASE0_FINDINGS_DECLARED_BEFORE_FREEZE": {
   "reconciliation": {"chromosomal_occurrences": 35140, "blocks": 21955, "mge_features": 32364,
                      "occurrence_window_positive": 16303, "block_positive": 6617,
                      "wrapped_circular_blocks": 57, "occurrences_in_wrapped_blocks": 107,
                      "truncated_blocks": 5, "unannotated_occurrences": 0,
                      "bioprojects_over_chromosomal_occurrences": 2248},
   "distance_convention_validated": {
     "test": "recompute the nearest qualifying MGE distance from raw coordinates and compare to "
             "the frozen arg_mge_neighbourhood.nearest_mge_distance_bp column",
     "result": "16,303 reconstructed positive and 18,837 censored, matching the frozen column "
               "exactly; 16,303 exact distance agreements; 0 mismatches (0.0000 percent)",
     "consequence": "the convention is confirmed: overlap is 0, otherwise the gap in bp between "
                    "nearest boundaries, scoped to the occurrence's own +/-10 kb window"},
   "CRITICAL_STRUCTURAL_FINDING": {
     "finding": "a BLOCK is a merged shared-context region that can be far larger than any single "
                "ARG's +/-10 kb window (block spans up to 62,349 bp were observed). A feature "
                "elsewhere in the same block can lie beyond 10 kb from a given ARG.",
     "evidence": "111 occurrences sit in an MGE-positive block yet have zero qualifying features "
                 "in their own window; their nearest in-block feature lies between 10,231 and "
                 "16,931 bp away. 16,414 block-positive minus 111 equals 16,303 window-positive, "
                 "which equals frozen class B exactly.",
     "consequence_for_this_protocol": "distances are computed STRICTLY within each occurrence's "
                                      "own +/-10 kb window. Those 111 occurrences are "
                                      "RIGHT-CENSORED at 10,000 bp even though an in-block "
                                      "feature exists beyond the window. Using it would "
                                      "manufacture a distance the frozen design declares "
                                      "unavailable and would breach the +/-20 kb NOT EVALUABLE "
                                      "rule.",
     "why_declared_here": "this was discovered during read-only feasibility, before any "
                          "group-specific outcome was computed, and is recorded so that the "
                          "censoring rule cannot be mistaken for a post-hoc choice"}},

 "population": {
   "occurrences": "all 35,140 PRIMARY chromosomal acquired-ARG occurrences (classes A and B)",
   "blocks": "all 21,955 frozen shared-context blocks",
   "exclusions": "none based on marker status",
   "sequence_availability": "windows provide sequence only through +/-10 kb",
   "censoring": "absence of a qualifying marker within the window is RIGHT-CENSORED at 10,000 bp "
                "and is never assigned an exact distance beyond 10 kb",
   "twenty_kb": "NOT EVALUABLE; no +/-20 kb quantity may be reported"},

 "distance_definition": {
   "overlap": "distance = 0 when the ARG interval and the MGE-feature interval overlap",
   "otherwise": "the number of base pairs between the nearest interval boundaries",
   "circular": "topology-aware coordinates; separation is min over the linear gap and the gap "
               "under +/- one replicon length, so the artificial linearised origin is never "
               "treated as a biological boundary",
   "truncated": "truncated-window flags are preserved and carried into every output row",
   "scope": "the nearest qualifying MGE feature within the available +/-10 kb sequence"},

 "endpoints": {
   "primary": "nearest qualifying MGE marker of EITHER type (IS/transposase evidence or "
              "integrase/integron evidence)",
   "secondary": ["nearest IS/transposase", "nearest integrase/integron",
                 "overlap at distance 0", "marker detection within 1, 2, 5 and 10 kb"],
   "marker_evidence_not_redefined": True,
   "feature_classes_as_frozen": {"IS_or_transposase": 29331, "integrase_or_integron": 3033,
                                 "total": 32364},
   "integrase_arm_is_secondary_because": "its evidence arm is sparse (3,033 of 32,364 features) "
                                         "and NM-V1 did not independently validate oriT or every "
                                         "NON-EVALUABLE boundary"},

 "primary_unit_and_weighting": {
   "estimand": "block-balanced occurrence distribution",
   "rule": "each shared block contributes total weight 1; if block b contains m ARG occurrences, "
           "each of those occurrences receives weight 1/m",
   "rationale": "a block containing several ARGs cannot dominate through multiplicity alone",
   "retention": "the distance of EVERY ARG occurrence is retained; a block is never replaced by "
                "its minimum distance in the primary estimand"},

 "sensitivity_estimands": {
   "S1": "occurrence-weighted; every ARG occurrence has weight 1",
   "S2": "one deterministic ARG occurrence per block, selected by the lowest SHA-256 of "
         "'<block_id>|<occurrence_id>|20260822'",
   "S3": "minimum ARG-to-MGE distance per block, labelled 'any-ARG block distance', which "
         "FAVOURS multi-ARG blocks and is reported as such",
   "S4": "excluding the five truncated blocks",
   "S5": "excluding circular-wrapped blocks",
   "S6": "IS/transposase-only endpoint"},

 "groups_and_contrasts": {
   "primary_groups": {
     "Acinetobacter baumannii": "organism_harmonized == 'Acinetobacter baumannii'",
     "Pseudomonas aeruginosa": "organism_harmonized == 'Pseudomonas aeruginosa'",
     "Klebsiella group": "genus == 'Klebsiella', exactly as defined in NM-V4 and NM-V4C"},
   "primary_contrasts": {"P1": "A. baumannii versus Klebsiella",
                          "P2": "A. baumannii versus P. aeruginosa"},
   "secondary_contrast": {"P3": "P. aeruginosa versus Klebsiella"},
   "all_other_species": "DESCRIPTIVE ONLY unless a numbered amendment is frozen before their "
                        "outcomes are examined"},

 "summary_measures": {
   "cumulative_incidence": "weighted F(d) = P(nearest MGE distance <= d) at 1,000, 2,000, 5,000 "
                           "and 10,000 bp",
   "full_curve": "weighted distance curve from 0 to 10,000 bp",
   "restricted_mean": "restricted mean distance-to-marker through 10 kb, treating "
                      "marker-negative occurrences as right-censored at 10 kb",
   "median": "reported only if estimable before 10 kb, otherwise 'not reached'",
   "auc": "area under the cumulative detection curve through 10 kb, as a descriptive "
          "spatial-proximity summary",
   "prohibited": "no ordinary mean after assigning marker-negative observations a fabricated "
                 "distance"},

 "uncertainty": {
   "primary": "BioProject-cluster bootstrap",
   "resamples": 2000, "seed": SEED,
   "procedure": "resample BioProjects with replacement; preserve every block and occurrence "
                "belonging to each sampled BioProject; RECALCULATE block weights inside every "
                "resample",
   "intervals": "percentile 95 percent CIs for curves and contrasts",
   "also_reported": "raw BioProject counts and effective 1/HHI",
   "prohibited": "occurrence-level binomial intervals as the primary uncertainty"},

 "multiplicity": {
   "primary_contrasts_frozen": 2, "primary_landmarks_frozen": [1000, 2000, 5000, 10000],
   "correction": "Holm across the two primary species contrasts, applied separately at each "
                 "prespecified landmark",
   "no_threshold_search": "no search for an optimal distance threshold is permitted",
   "hierarchy": "curves and effect sizes are primary; p-values are secondary"},

 "required_interpretation": {
   "A": "does the A. baumannii result appear across the distance curve rather than only at 10 kb",
   "B": "is it concentrated at very short distances, consistent with tight ARG-MGE architecture",
   "C": "does it remain after block balancing and BioProject resampling",
   "D": "does it remain under the deterministic one-ARG-per-block sensitivity",
   "E": "is it predominantly driven by IS/transposase rather than the smaller integrase component",
   "prohibited_claims": ["intact transposons", "demonstrated mobilization", "HGT", "transfer",
                          "phenotype"]},

 "stop_conditions": [
   "35,140 occurrences do not reconcile exactly",
   "21,955 blocks do not reconcile exactly",
   "any occurrence maps ambiguously to multiple blocks",
   "more than five truncated blocks are found",
   "circular wrapping cannot be reconstructed consistently",
   "frozen MGE calls or class assignments differ from their recorded hashes",
   "the primary block-balanced weights do not sum exactly to the number of represented blocks",
   "more than 1 percent of exact distances cannot be reconstructed",
   "a post-outcome change to endpoint, group or threshold becomes necessary"],
 "amendment_rule": "any required post-outcome change needs a numbered amendment; the protocol is "
                   "never silently modified",

 "governance": {"execution": "local only", "server_use": False, "installation_or_download": False,
                "dependencies_present": {"python": "3.14.6", "numpy": "2.5.0",
                                          "matplotlib": "3.11.1", "scipy": "1.18.0"},
                "wall_time_budget_hours": 8,
                "prohibited": ["Paper 2", "PlasmidCall", "MicContext", "commit", "push",
                               "modification of frozen artefacts", "rewriting the manuscript"]},
 "input_digests": dig,
    }
    json.dump(Q, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("%s\n  wrote %s" % (VERSION, a.out))
    print("  FROZEN PROTOCOL SHA-256: %s" % sha(a.out))


if __name__ == "__main__":
    main()
