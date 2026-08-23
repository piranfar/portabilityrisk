"""NM-V1 ground-truth construction rules -- frozen BEFORE aggregate performance is computed.

The frozen adjudicated labels, not any single tool output, are the NM-V1 biological ground
truth. This file fixes which cases may receive an automatic label without human review, which
must go to blinded adjudication, and how the quality-control arm is drawn -- all before any
sensitivity, specificity or agreement statistic is calculated.

The QC arm exists because concordance between two tools is not truth: both share HMMER,
Prodigal and a protein-homology view of transposases, so they can agree and both be wrong.
Sampling concordant blocks into the blinded set is the only way to detect that.
"""
import argparse, datetime, hashlib, json, os, sys

VERSION = "nmv1_freeze_groundtruth_rules_v1.0.0"
DESIGN_SHA = "c2aea6cb583c24b997ab376861acc600295e94d59bfa0ef2d55cf2bcc424bb20"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


R = {
 "rules": "NMV1_GROUNDTRUTH_CONSTRUCTION_RULES",
 "version": "1.0.0", "builder": VERSION,
 "frozen_before_aggregate_performance_was_computed": True,
 "nmv1_design_sha256": DESIGN_SHA,
 "adjudicator": {
   "name": "Vahhab Piranfar",
   "eligibility_basis": "did not develop the original MGE annotation pipeline or the NM-V1 "
                        "scoring code, and has not inspected block-level validation outcomes",
   "authorised": "owner, D-NM3, 2026-08-21"},

 "principle": "The frozen adjudicated labels are the ground truth. No single tool output is "
              "truth. ISEScan and IntegronFinder are INDEPENDENT STRUCTURAL EVIDENCE, not "
              "gold standards: both share HMMER and Prodigal with the method under test.",

 "evidence_fields_per_block": [
   "original HMM annotation (IS and integron counts)",
   "ISEScan annotation: element count, family, cluster, complete vs partial type",
   "IntegronFinder annotation: complete integron, CALIN, In0 counts",
   "genomic coordinates of the block on its replicon",
   "element boundaries where reported (isBegin, isEnd, isLen)",
   "inverted-repeat evidence (irId, irLen, tir sequences)",
   "transposase evidence (orfBegin, orfEnd, orfLen, E-value)",
   "integron structural evidence (attC sites, integrase)",
   "topology and truncation status (circular, linear, wrapped, truncated)",
   "tool execution status (completed, failed, timeout, no-output)"],

 "automatic_label_rules": {
   "applies_only_to": "concordant high-confidence structural findings",
   "AUTO_IS_MOBILE": {
     "condition": "HMM IS-positive AND ISEScan reports at least one element with type = c "
                  "(complete: both terminal inverted repeats resolved and a full transposase ORF)",
     "label": "IS-associated supported",
     "maps_to": "chromosomal_mobile",
     "why_automatic": "a complete IS with resolved TIRs detected by a boundary-based method, "
                      "agreeing with protein homology, is the least ambiguous evidence available"},
   "AUTO_INTEGRON": {
     "condition": "HMM integron-positive AND IntegronFinder reports at least one COMPLETE "
                  "integron (integrase plus at least one attC site)",
     "label": "integron-associated supported",
     "maps_to": "chromosomal_mobile"},
   "AUTO_QUIESCENT": {
     "condition": "HMM negative for both marker classes AND ISEScan reports zero elements of "
                  "any type AND IntegronFinder reports zero complete, zero CALIN and zero In0",
     "label": "chromosomal-quiescent supported",
     "maps_to": "chromosomal_quiescent"},
   "no_other_automatic_label_exists": True},

 "mandatory_adjudication": {
   "rule": "every block below enters blinded biological adjudication; none may be auto-labelled",
   "categories": {
     "DISCORDANT_IS": "HMM IS status disagrees with ISEScan positivity",
     "DISCORDANT_INTEGRON": "HMM integron status disagrees with IntegronFinder positivity",
     "PARTIAL_ONLY": "ISEScan reports elements but none with type = c; partial elements alone "
                     "do not establish a mobile context",
     "INCOMPLETE_INTEGRON": "IntegronFinder reports CALIN or In0 but no complete integron",
     "BOUNDARY_AMBIGUOUS": "block is truncated or wrapped_circular, so an element may extend "
                           "beyond the retrieved sequence",
     "TOOL_FAILURE": "any tool failed, timed out, or produced no parsable output",
     "REFERENCE_DISAGREEMENT": "ISEScan and IntegronFinder imply different architectures for "
                               "the same block"}},

 "quality_control_arm": {
   "purpose": "test whether agreement between the two reference tools conceals a shared error. "
              "Both use HMMER and Prodigal, so concordance is not independence.",
   "rule": "a prespecified random sample of blocks that WOULD have been auto-labelled is added "
           "to the blinded adjudication set, indistinguishable from the discordant cases",
   "allocation": {"AUTO_IS_MOBILE": 40, "AUTO_QUIESCENT": 40, "AUTO_INTEGRON": 20},
   "seed": 20260821,
   "blinding": "QC blocks carry the same opaque token format and are shuffled together with "
               "the discordant cases; the adjudicator cannot tell them apart",
   "interpretation": "if adjudication overturns auto-labels at a materially higher rate than "
                     "chance, the automatic rules are unsafe and every auto-label must be "
                     "reconsidered before any performance figure is reported"},

 "permitted_adjudication_outcomes": [
   "chromosomal_mobile_supported", "chromosomal_quiescent_supported",
   "integron_associated_supported", "IS_associated_supported",
   "multiple_MGE_evidence_supported", "neither_classification_supported",
   "biologically_indeterminate"],

 "adjudication_rubric": {
   "instruction": "Judge only the structural evidence shown. Do not infer from gene identity, "
                  "organism or study; none is provided.",
   "IS_associated_supported": "a transposase ORF with at least one resolved terminal inverted "
                              "repeat, or a complete element reported by a boundary-based method",
   "integron_associated_supported": "an integrase with at least one attC site, or a complete "
                                    "integron structure",
   "multiple_MGE_evidence_supported": "both IS and integron evidence present in the block",
   "chromosomal_mobile_supported": "any of the three above; the block sits in mobile context",
   "chromosomal_quiescent_supported": "no credible mobile-element evidence within the block",
   "neither_classification_supported": "evidence is present but contradicts both readings",
   "biologically_indeterminate": "the evidence shown cannot decide; USE THIS RATHER THAN "
                                 "GUESSING. An honest indeterminate is more useful than a "
                                 "forced call and is reported as its own category.",
   "boundary_cases": "if an element appears truncated at a block edge, judge only what is "
                     "visible and mark indeterminate if the visible part is insufficient"},

 "fields_withheld_from_the_adjudicator": [
   "species", "BioProject", "assembly accession", "ARG identity", "gene family",
   "original portability class", "sampling stratum", "tool names (presented as Method X, "
   "Method Y, Method Z)", "any aggregate result"],
 "fields_shown_because_essential_for_biological_interpretation": [
   "sequence segment", "topology", "coordinates within the block", "element boundaries",
   "inverted-repeat sequences and identity", "ORF coordinates and length", "E-values",
   "attC and integrase evidence", "truncation status"],

 "unblinding_discipline": {
   "key_written_and_hashed_before_delivery": True,
   "key_not_exposed_until": "the completed adjudication file has been returned, frozen and hashed",
   "original_decisions_preserved": "adjudicator decisions are never rewritten after unblinding; "
                                   "a rule correction requires a dated amendment and the "
                                   "original is retained"},

 "performance_reporting_rules": {
   "ground_truth": "frozen adjudicated labels, plus auto-labels that survive the QC arm",
   "design_weighting": "mandatory; strata were sampled at 2.6 to 100 per cent, so unweighted "
                       "estimates are biased and are reported only as secondary descriptives",
   "confidence_intervals": "stratified bootstrap, 2000 resamples, seed 20260821",
   "arms_reported_separately_before_any_combined_class_B_figure": True,
   "prohibited": ["calling homology-only detections false positives outside the adjudicated set",
                  "pooling the IS and integron arms before reporting them separately",
                  "overwriting the frozen 74,349-occurrence dataset - any recomputation is an "
                  "additive validation layer with its own name"]},

 "distinctions_that_must_be_preserved_in_reporting": {
   "validated_genomic_portability_architecture": "what NM-V1 can establish",
   "predicted_mobility_potential": "what MOB-suite and marker presence indicate",
   "experimentally_observed_horizontal_transfer": "NOT measured anywhere in this programme",
   "note": "the absence of a conjugation experiment does not invalidate a correctly validated "
           "genomic architecture and must not be presented as if it does"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if sha256_file(a.design) != DESIGN_SHA:
        print("REFUSING: NM-V1 design digest mismatch"); sys.exit(1)
    if os.path.exists(a.out):
        print("REFUSING: %s exists; frozen rules are never overwritten" % a.out); sys.exit(1)
    R["frozen_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    json.dump(R, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("%s\n  wrote %s" % (VERSION, a.out))
    print("  FROZEN RULES SHA-256: %s" % sha256_file(a.out))
    print("  adjudicator: %s" % R["adjudicator"]["name"])
    print("  auto-label rules: %d | mandatory adjudication categories: %d | QC arm: %d blocks"
          % (3, len(R["mandatory_adjudication"]["categories"]),
             sum(R["quality_control_arm"]["allocation"].values())))


if __name__ == "__main__":
    main()
